"""
``labw_utils.commonutils.str_utils`` -- String utilities
    
.. versionadded:: 1.0.2
"""

from __future__ import annotations

import string

from labw_utils.typing_importer import Any, Optional, Dict, Mapping


def to_dict(
    in_str: str,
    field_sep: str = "\t",
    record_sep: str = "\n",
    quotation_mark: Optional[str] = None,
    resolve_str: bool = True,
) -> Mapping[str, Any]:
    """
    A simple parser to get key-value pairs to a dictionary.

    Key: string. Value: String, float or int.

    >>> input_str = '\\nCPU:\\t2\\nMEM:\\t5.1\\nPCIE:\\t3rd Gen\\nGRAPHICS:\\t"UHD630\\tRTX2070"\\nUSB: "3.1"\\nOthers:::info'
    >>> to_dict(input_str, field_sep=':', record_sep='\\n', quotation_mark="\\'\\"", resolve_str=True)
    {'CPU': 2, 'MEM': 5.1, 'PCIE': '3rd Gen', 'GRAPHICS': 'UHD630\\tRTX2070', 'USB': 3.1, 'Others': 'info'}
    >>> to_dict(input_str, field_sep=':', record_sep='\\n', quotation_mark="\\'\\"", resolve_str=False)
    {'CPU': '2', 'MEM': '5.1', 'PCIE': '3rd Gen', 'GRAPHICS': 'UHD630\\tRTX2070', 'USB': '3.1', 'Others': 'info'}
    >>> to_dict(input_str, field_sep=':', record_sep='\\n', quotation_mark=None, resolve_str=False)
    {'CPU': '2', 'MEM': '5.1', 'PCIE': '3rd Gen', 'GRAPHICS': '"UHD630\\tRTX2070"', 'USB': '"3.1"', 'Others': 'info'}
    >>> to_dict(input_str, field_sep=':', record_sep='\\n', quotation_mark=None, resolve_str=True)
    {'CPU': 2, 'MEM': 5.1, 'PCIE': '3rd Gen', 'GRAPHICS': '"UHD630\\tRTX2070"', 'USB': '"3.1"', 'Others': 'info'}
    >>> to_dict(
    ... 'CPU:\\t2\\nMEM:\\t5.1\\nGRAPHICS:\\tUHD630\\nGRAPHICS:\\tRTX2070',
    ... field_sep=':', record_sep='\\n', quotation_mark=None, resolve_str=False
    ... )
    {'CPU': '2', 'MEM': '5.1', 'GRAPHICS': ['UHD630', 'RTX2070']}

    :param in_str: Input string
    :param field_sep: Field separator, the FS variable in AWK programming language.
    :param record_sep: Record separator, the RS variable in AWK programming language.
    :param quotation_mark: If the key and value is not quoted, pass ``None``.
                           If quoted by single quote, pass ``'\''``.
                           If quoted by double quote, pass ``'\"'``.
                           If quoted by single and double quote, pass ``'\'\"'``.
                           Will not parse quoted by triple quotes.
    :param resolve_str: Whether to resolve strings to float or int.

    .. versionadded:: 1.0.2
    """
    retd: Dict[str, Any] = {}
    record_val: Any
    in_str_by_record = in_str.split(record_sep)
    for record in in_str_by_record:
        record = record.strip(string.whitespace + field_sep)
        lr = len(record)
        first_field_sep_pos = record.find(field_sep)
        if first_field_sep_pos == -1:
            continue
        record_key = record[0:first_field_sep_pos].rstrip()
        while first_field_sep_pos != lr:
            if record[first_field_sep_pos] == field_sep:
                first_field_sep_pos += 1
            else:
                record_val = record[first_field_sep_pos:].lstrip()
                break
        else:
            record_val = ""
        if quotation_mark is not None:
            record_key = record_key.strip(quotation_mark)
            record_val = record_val.strip(quotation_mark)
        if resolve_str:
            try:
                if "." in record_val:
                    record_val = float(record_val)
                else:
                    record_val = int(record_val)
            except ValueError:
                pass
        if record_key in retd:
            if isinstance(retd[record_key], list):
                retd[record_key].append(record_val)
            else:
                retd[record_key] = [retd[record_key], record_val]
        else:
            retd[record_key] = record_val
    return retd
