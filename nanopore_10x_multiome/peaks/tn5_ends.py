"""
claude-sonnet-4.6
"""

import pysam


# ---------------------------------------------------------------------------
# CIGAR operation classification
#
# SAM CIGAR operations are numbered 0–8:
#   0 M  alignment match (can be either a sequence match or mismatch)
#   1 I  insertion relative to the reference
#   2 D  deletion from the reference
#   3 N  skipped region in the reference (intron in RNA; treated like D here)
#   4 S  soft clip (bases present in query, not aligned)
#   5 H  hard clip (bases absent from query sequence record)
#   6 P  padding (silent deletion from padded reference)
#   7 =  sequence match
#   8 X  sequence mismatch
#
# "Query-consuming" ops advance the position in the read sequence.
# "Reference-consuming" ops advance the position on the reference genome.
# An op can be both (M, =, X), one, or neither (H, P).
# ---------------------------------------------------------------------------

_QUERY_CONSUMING = {0, 1, 4, 7, 8}   # M, I, S, =, X
_REF_CONSUMING = {0, 2, 3, 7, 8}     # M, D, N, =, X


def _trim_cigar_from_left(cigar_tuples, n):
    """Trim a CIGAR string to retain the first *n* query bases.

    Used to build the left-hand (5'-end) output read, which starts at the
    same reference position as the original read.

    Deletion (D) and skip (N) operations that appear between query-consuming
    ops are buffered and only committed once a following query-consuming op
    confirms they are interior to the kept region.  Any D/N ops that would
    appear at the *trailing* edge of the trimmed CIGAR are discarded; ending
    on a deletion would produce an invalid alignment.

    Parameters
    ----------
    cigar_tuples : list of (int, int)
        CIGAR as returned by ``pysam.AlignedSegment.cigartuples``, i.e. a list
        of ``(operation_code, length)`` pairs in left-to-right genome order.
    n : int
        Number of query bases to keep from the left (5') end.

    Returns
    -------
    new_cigar : list of (int, int)
        Trimmed CIGAR tuples covering exactly *n* query bases (or fewer if the
        read is shorter than *n*).
    ref_consumed : int
        Total reference bases spanned by *new_cigar*.  Adding this to the
        original ``reference_start`` gives the reference end of the output read.
    """
    new_cigar = []
    remaining = n        # query bases still to be assigned to the output CIGAR
    ref_consumed = 0
    pending = []         # D/N ops buffered until we know if they are interior

    for op, length in cigar_tuples:
        if remaining == 0:
            break

        if op in _QUERY_CONSUMING:
            # Flush any preceding D/N ops — they are interior to the kept region.
            for pop, plen in pending:
                new_cigar.append((pop, plen))
                ref_consumed += plen
            pending = []

            take = min(length, remaining)
            new_cigar.append((op, take))
            remaining -= take
            if op in _REF_CONSUMING:
                ref_consumed += take

        elif op in _REF_CONSUMING:
            # D or N: buffer it; it is only interior if a query-consuming op follows.
            pending.append((op, length))

    # Any ops remaining in `pending` would trail the last kept query base —
    # drop them so the CIGAR does not end on a deletion.
    return new_cigar, ref_consumed


def _trim_cigar_from_right(cigar_tuples, n):
    """Trim a CIGAR string to retain the last *n* query bases.

    Used to build the right-hand (3' fragment end / 5' cut-site) output read.
    The trimmed CIGAR must be re-anchored: the caller computes the new
    ``reference_start`` by subtracting *ref_consumed* from the original read's
    ``reference_end``.

    The algorithm mirrors :func:`_trim_cigar_from_left` but iterates the CIGAR
    in reverse.  D/N ops are buffered while travelling right-to-left, then
    inserted *after* (i.e., to the right of, in genome order) the preceding
    query-consuming op once that op is committed.  Trailing D/N ops at the
    leftmost edge are discarded.

    Parameters
    ----------
    cigar_tuples : list of (int, int)
        CIGAR as returned by ``pysam.AlignedSegment.cigartuples``.
    n : int
        Number of query bases to keep from the right (3') end.

    Returns
    -------
    new_cigar : list of (int, int)
        Trimmed CIGAR tuples covering exactly *n* query bases (or fewer).
    ref_consumed : int
        Total reference bases spanned by *new_cigar*.
    """
    new_cigar = []
    remaining = n        # query bases still to be assigned
    ref_consumed = 0
    pending = []         # D/N ops buffered (accumulated right-to-left)

    for op, length in reversed(cigar_tuples):
        if remaining == 0:
            break

        if op in _QUERY_CONSUMING:
            take = min(length, remaining)
            # Prepend so the final list is in left-to-right genome order.
            new_cigar.insert(0, (op, take))
            remaining -= take
            if op in _REF_CONSUMING:
                ref_consumed += take

            # Flush buffered D/N ops to the RIGHT of the current op (positions
            # 1, 2, … in new_cigar) because they were encountered before it
            # in reverse traversal but are genomically downstream of it.
            for i, (pop, plen) in enumerate(pending, 1):
                new_cigar.insert(i, (pop, plen))
                ref_consumed += plen
            pending = []

        elif op in _REF_CONSUMING:
            # D or N encountered before any query-consuming op in this reversed
            # pass — buffer it; it may be interior if a query op follows.
            pending.insert(0, (op, length))

    # Remaining `pending` entries would be at the leftmost edge — discard them.
    return new_cigar, ref_consumed


def _make_single_end_read(read, n_bp, header):
    """Produce two trimmed single-end reads from one paired-end alignment.

    Each output read carries the *n_bp* query bases nearest to a Tn5 cut site:

    * ``outs[0]`` — left-hand read: first *n_bp* bases, anchored at the read's
      original ``reference_start`` (the 5' cut site for forward-strand reads).
    * ``outs[1]`` — right-hand read: last *n_bp* bases, anchored so that its
      right edge aligns with the read's original ``reference_end`` (the 5' cut
      site for reverse-strand reads).

    Paired-end SAM flags (``0x1 paired``, ``0x2 proper pair``, ``0x8 mate
    unmapped``, ``0x20 mate reverse``, ``0x40 read1``, ``0x80 read2``) are
    cleared so downstream tools treat the output reads as single-end.  Mate
    CIGAR (MC), mate score (ms), and mate mapping quality (MQ) tags are also
    dropped as they no longer apply.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A single mapped paired-end read with a valid CIGAR string.
    n_bp : int
        Number of bases to retain from each end of the read sequence.
    header : pysam.AlignmentHeader
        BAM header to attach to the output segments.

    Returns
    -------
    tuple of (pysam.AlignedSegment, pysam.AlignedSegment)
        ``(left_read, right_read)`` as described above.
    None
        If the read has no sequence or is shorter than ``2 * n_bp`` (not
        enough bases to produce two non-overlapping trimmed reads).
    """
    # Require at least 2×n_bp bases so the two trimmed reads do not overlap.
    if read.query_sequence is None or len(read.query_sequence) < (n_bp * 2):
        return None

    # Allocate two output segments sharing the same header.
    # outs[0] = left-hand (5' of read), outs[1] = right-hand (3' of read).
    outs = (
        pysam.AlignedSegment(header),
        pysam.AlignedSegment(header)
    )

    for out in outs:
        out.query_name = read.query_name
        out.mapping_quality = read.mapping_quality
        out.reference_id = read.reference_id

        # Clear mate information — these are now unpaired single-end reads.
        out.next_reference_id = -1
        out.next_reference_start = -1
        out.template_length = 0

        # Drop mate-specific tags that are meaningless for single-end reads.
        out.tags = [(t, v) for t, v in read.tags if t not in {'MC', 'ms', 'MQ'}]

        # Clear paired-end flag bits while preserving strand, secondary,
        # supplementary, QC-fail, and duplicate flags.
        flag = read.flag
        for bit in (0x1,   # read paired
                    0x2,   # read mapped in proper pair
                    0x8,   # mate unmapped
                    0x20,  # mate reverse strand
                    0x40,  # read is read 1
                    0x80): # read is read 2
            flag &= ~bit
        out.flag = flag

    # Assign trimmed sequences and base qualities.
    outs[0].query_sequence = read.query_sequence[:n_bp]
    outs[1].query_sequence = read.query_sequence[-n_bp:]

    if read.query_qualities is not None:
        outs[0].query_qualities = read.query_qualities[:n_bp]
        outs[1].query_qualities = read.query_qualities[-n_bp:]

    # Build trimmed CIGAR strings for each output read.
    outs[0].cigar, _ = _trim_cigar_from_left(read.cigartuples, n_bp)

    # _trim_cigar_from_right returns the number of reference bases consumed by
    # the trimmed CIGAR; subtract from reference_end to get the correct start.
    outs[1].cigar, ref_consumed = _trim_cigar_from_right(read.cigartuples, n_bp)

    # Left read starts at the same reference position as the original read.
    outs[0].reference_start = read.reference_start
    # Right read ends at the same reference position as the original read.
    outs[1].reference_start = read.reference_end - ref_consumed

    return outs


def paired_end_to_single_end_trim(input_bam, output_bam, n_bp=25, pbar=False):
    """Convert an aligned paired-end BAM to a 5'-trimmed single-end BAM.

    For each mapped read in *input_bam* the function emits two single-end
    output reads (one per Tn5 cut site) that each retain only the *n_bp* bases
    nearest to the cut site.  This representation is the standard input for
    ATAC-seq peak callers such as MACS2/3, which use the read 5' ends as
    proxies for Tn5 insertion positions.

    Reads that are unmapped, lack a CIGAR string, or are shorter than
    ``2 * n_bp`` bases are silently dropped.  The output BAM inherits the
    header (including sequence dictionary and read groups) from *input_bam*.

    Parameters
    ----------
    input_bam : str or path-like
        Path to the aligned paired-end BAM file.  Must be sorted and the index
        is not required; the file is read sequentially.
    output_bam : str or path-like
        Path for the output single-end BAM.  The file is created (or
        overwritten) and written in BAM format.
    n_bp : int, optional
        Number of bases to retain from the 5' end of each read.  Defaults to
        25, which is sufficient for most peak callers while keeping file size
        small.
    pbar : bool, optional
        If ``True``, display a :mod:`tqdm` progress bar reporting the number of
        reads processed.  Requires :mod:`tqdm` to be installed.  Defaults to
        ``False``.

    Examples
    --------
    >>> paired_end_to_single_end_trim("aligned.bam", "trimmed.bam", n_bp=25)

    >>> # With a progress bar:
    >>> paired_end_to_single_end_trim("aligned.bam", "trimmed.bam", pbar=True)
    """
    with pysam.AlignmentFile(input_bam, 'rb') as bam_in:
        with pysam.AlignmentFile(output_bam, 'wb', template=bam_in) as bam_out:

            if pbar:
                import tqdm
                reads = tqdm.tqdm(bam_in, total=bam_in.mapped + bam_in.unmapped)
            else:
                reads = bam_in

            for read in reads:
                # Skip reads without a valid alignment — they carry no cut-site
                # information and cannot be CIGAR-trimmed.
                if read.is_unmapped or read.cigartuples is None:
                    continue

                outs = _make_single_end_read(read, n_bp, bam_in.header)
                if outs is not None:
                    for out in outs:
                        bam_out.write(out)
