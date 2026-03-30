"""
Property-based tests for quant.tsv Name Normalization.

# Feature: neoantigen-workflow, Property 9: quant.tsv Name Normalization

Validates: Requirements 4a.2, 4a.3, 4a.4

For any quant.sf file, the quant.tsv produced by QUANT_TSV_NORMALIZE must
satisfy:
  (a) every value in the Name column contains no '|' character, and
  (b) all other columns (Length, EffectiveLength, TPM, NumReads) are
      byte-for-byte identical to the corresponding values in quant.sf.

The normalization logic implemented in modules/local/quant_tsv_normalize/main.nf:
    parts[0] = parts[0].split('|')[0]
"""

from hypothesis import given, settings
from hypothesis import strategies as st

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

QUANT_SF_HEADER = "Name\tLength\tEffectiveLength\tTPM\tNumReads"

# ---------------------------------------------------------------------------
# Normalization logic (mirrors the Nextflow inline Python)
# ---------------------------------------------------------------------------


def normalize_name(name: str) -> str:
    """Strip everything after the first '|' in a transcript name."""
    return name.split("|")[0]


def normalize_quant_sf(content: str) -> str:
    """
    Apply QUANT_TSV_NORMALIZE logic to quant.sf content.

    Reads lines, strips '|'-suffix from Name column (skipping header),
    returns the resulting quant.tsv content as a string.
    """
    lines = content.splitlines()
    out_lines = []
    for line in lines:
        parts = line.split("\t")
        if parts[0] != "Name":  # skip header
            parts[0] = normalize_name(parts[0])
        out_lines.append("\t".join(parts))
    return "\n".join(out_lines)


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Plain transcript ID (no '|'), e.g. ENST00000456328.2
plain_name_st = st.text(
    alphabet=st.characters(
        whitelist_categories=("Lu", "Ll", "Nd"),
        whitelist_characters="._-",
    ),
    min_size=1,
    max_size=30,
)

# Gencode pipe-delimited name, e.g. ENST00000456328.2|ENSG00000223972.5|...
gencode_name_st = st.builds(
    lambda parts: "|".join(parts),
    st.lists(
        st.text(
            alphabet=st.characters(
                whitelist_categories=("Lu", "Ll", "Nd"),
                whitelist_characters="._-",
            ),
            min_size=1,
            max_size=20,
        ),
        min_size=2,
        max_size=5,
    ),
)

# Mix of plain and pipe-delimited names
transcript_name_st = st.one_of(plain_name_st, gencode_name_st)

# A single quant.sf data line (5 tab-separated columns)
data_line_st = st.builds(
    lambda name, length, eff_length, tpm, num_reads: (
        f"{name}\t{length}\t{eff_length:.6f}\t{tpm:.6f}\t{num_reads:.6f}"
    ),
    name=transcript_name_st,
    length=st.integers(min_value=100, max_value=10000),
    eff_length=st.floats(min_value=1.0, max_value=10000.0, allow_nan=False, allow_infinity=False),
    tpm=st.floats(min_value=0.0, max_value=1e6, allow_nan=False, allow_infinity=False),
    num_reads=st.floats(min_value=0.0, max_value=1e8, allow_nan=False, allow_infinity=False),
)

# Complete quant.sf content: header + 1–100 data lines
quant_sf_content_st = st.builds(
    lambda lines: QUANT_SF_HEADER + "\n" + "\n".join(lines),
    lines=st.lists(data_line_st, min_size=1, max_size=100),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def parse_quant(content: str) -> tuple[str, list[list[str]]]:
    """
    Parse quant file content into (header_line, list_of_column_lists).
    Returns the raw header string and each data line split by tab.
    """
    lines = content.splitlines()
    header = lines[0] if lines else ""
    data_rows = [line.split("\t") for line in lines[1:] if line.strip()]
    return header, data_rows


# ---------------------------------------------------------------------------
# Property 9a: No '|' in Name column after normalization
# ---------------------------------------------------------------------------


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p9_no_pipe_in_name_column_after_normalization(content):
    """
    **Validates: Requirements 4a.2, 4a.3**

    Property 9: quant.tsv Name Normalization.
    After normalization, no value in the Name column of quant.tsv may
    contain a '|' character.

    # Feature: neoantigen-workflow, Property 9: quant.tsv Name Normalization
    """
    normalized = normalize_quant_sf(content)
    _, data_rows = parse_quant(normalized)

    assert len(data_rows) >= 1, "quant.tsv must have at least one data row"

    for i, row in enumerate(data_rows):
        name = row[0]
        assert "|" not in name, (
            f"Row {i + 1}: Name column still contains '|' after normalization. "
            f"Got: {name!r}"
        )


# ---------------------------------------------------------------------------
# Property 9b: Other columns are byte-for-byte identical to input
# ---------------------------------------------------------------------------


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p9_other_columns_unchanged_after_normalization(content):
    """
    **Validates: Requirements 4a.4**

    Property 9: quant.tsv Name Normalization.
    All columns other than Name (Length, EffectiveLength, TPM, NumReads)
    must be byte-for-byte identical between quant.sf and quant.tsv.

    # Feature: neoantigen-workflow, Property 9: quant.tsv Name Normalization
    """
    normalized = normalize_quant_sf(content)

    _, input_rows = parse_quant(content)
    _, output_rows = parse_quant(normalized)

    assert len(input_rows) == len(output_rows), (
        f"Row count changed: input has {len(input_rows)} rows, "
        f"output has {len(output_rows)} rows."
    )

    for i, (in_row, out_row) in enumerate(zip(input_rows, output_rows)):
        # Columns 1–4: Length, EffectiveLength, TPM, NumReads
        for col_idx, col_name in enumerate(
            ["Length", "EffectiveLength", "TPM", "NumReads"], start=1
        ):
            assert in_row[col_idx] == out_row[col_idx], (
                f"Row {i + 1}, column '{col_name}' (index {col_idx}) changed "
                f"after normalization. "
                f"Input: {in_row[col_idx]!r}, Output: {out_row[col_idx]!r}"
            )


# ---------------------------------------------------------------------------
# Property 9c: Plain names (no '|') are left unchanged
# ---------------------------------------------------------------------------


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p9_plain_names_unchanged_after_normalization(content):
    """
    **Validates: Requirements 4a.3**

    Property 9: quant.tsv Name Normalization.
    If a Name value contains no '|', it must be identical in quant.tsv.

    # Feature: neoantigen-workflow, Property 9: quant.tsv Name Normalization
    """
    normalized = normalize_quant_sf(content)

    _, input_rows = parse_quant(content)
    _, output_rows = parse_quant(normalized)

    for i, (in_row, out_row) in enumerate(zip(input_rows, output_rows)):
        original_name = in_row[0]
        normalized_name = out_row[0]
        if "|" not in original_name:
            assert normalized_name == original_name, (
                f"Row {i + 1}: Plain name (no '|') was modified. "
                f"Input: {original_name!r}, Output: {normalized_name!r}"
            )


# ---------------------------------------------------------------------------
# Property 9d: Pipe-delimited names are truncated to the first segment
# ---------------------------------------------------------------------------


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p9_pipe_delimited_names_truncated_to_first_segment(content):
    """
    **Validates: Requirements 4a.2**

    Property 9: quant.tsv Name Normalization.
    For any Name value containing '|', the normalized value must equal
    exactly the substring before the first '|'.

    # Feature: neoantigen-workflow, Property 9: quant.tsv Name Normalization
    """
    normalized = normalize_quant_sf(content)

    _, input_rows = parse_quant(content)
    _, output_rows = parse_quant(normalized)

    for i, (in_row, out_row) in enumerate(zip(input_rows, output_rows)):
        original_name = in_row[0]
        normalized_name = out_row[0]
        if "|" in original_name:
            expected = original_name.split("|")[0]
            assert normalized_name == expected, (
                f"Row {i + 1}: Pipe-delimited name not correctly truncated. "
                f"Input: {original_name!r}, Expected: {expected!r}, Got: {normalized_name!r}"
            )
