"""
``labw_utils.stdlib_helper.argparse_helper`` -- Argument parser with enhanced help formatter

Following is an example using "normal" formatter:

>>> import sys
>>> import pytest
>>> import doctest
>>> parser = argparse.ArgumentParser(prog="prog", description="description")
>>> _ = parser.add_argument("p", type=int, help="p-value")
>>> _ = parser.add_argument("-i", required=True, type=str, help="input filename")
>>> _ = parser.add_argument("-o", required=False, type=str, help="output filename", default="/dev/stdout")
>>> _ = parser.add_argument("--flag", action="store_true", help="flag")

Please notice that the default format differs between Python 3.9 and 3.10.

Below is an example on how it would show on Python <= 3.9:

>>> if sys.version_info > (3, 9):
...     pytest.skip()
... else:
...     print(parser.format_help())
usage: prog [-h] -i I [-o O] [--flag] p
<BLANKLINE>
description
<BLANKLINE>
positional arguments:
  p           p-value
<BLANKLINE>
optional arguments:
  -h, --help  show this help message and exit
  -i I        input filename
  -o O        output filename
  --flag      flag
<BLANKLINE>

Following is an example using enhanced formatter:

>>> parser = ArgumentParserWithEnhancedFormatHelp(prog="prog", description="description")
>>> _ = parser.add_argument("p", type=int, help="p-value")
>>> _ = parser.add_argument("-i", required=True, type=str, help="input filename")
>>> _ = parser.add_argument("-o", required=False, type=str, help="output filename", default="/dev/stdout")
>>> _ = parser.add_argument("--flag", action="store_true", help="flag")
>>> print(parser.format_help())
description
<BLANKLINE>
SYNOPSIS: prog [-h] -i I [-o O] [--flag] p
<BLANKLINE>
PARAMETERS:
  p
              [REQUIRED] Type: int; No defaults
              p-value
<BLANKLINE>
OPTIONS:
  -h, --help
              [OPTIONAL]
              show this help message and exit
  -i I
              [REQUIRED] Type: str; No defaults
              input filename
  -o O
              [OPTIONAL] Type: str; Default: /dev/stdout
              output filename
  --flag
              [OPTIONAL] Default: False
              flag
<BLANKLINE>

Following is an example on how it deals with :py:class:`enum.Enum`:

>>> class SampleEnum(enum.Enum):
...     \"\"\"SampleEnum docstring\"\"\"
...     A = 1
...     B = 2
...
...     @classmethod
...     def from_name(cls, in_name: str) -> \'SampleEnum\':
...         for choices in cls:
...             if choices.name == in_name:
...                 return choices
...         raise TypeError
>>> SampleEnum.A.__doc__ = "A doc"
>>> SampleEnum.B.__doc__ = None
>>> parser = ArgumentParserWithEnhancedFormatHelp(prog="prog", description="description")
>>> _ = parser.add_argument(
...     "-e",
...     required=False,
...     type=SampleEnum.from_name,
...     choices=SampleEnum,
...     help="output filename",
...     default=SampleEnum.A
... )
>>> print(parser.format_help())
description
<BLANKLINE>
SYNOPSIS: prog [-h] [-e {A,B}]
<BLANKLINE>
OPTIONS:
  -h, --help
              [OPTIONAL]
              show this help message and exit
  -e {A,B}
              [OPTIONAL] Type: SampleEnum; Default: A
              output filename
              CHOICES:
                  A -- A doc
                  B -- None
<BLANKLINE>

>>> parser.parse_args(["-e", "A"])
Namespace(e=<SampleEnum.A: 1>)

.. versionadded:: 1.0.2
.. versionchanged:: 1.0.3
    Support of enumns added.
    Setting both ``required`` and ``default`` would generate :py:class:`ValueError`.

.. warning:: Support of enumns highly experimental!

.. todo:: Markdown support.
"""

__all__ = ("ArgumentParserWithEnhancedFormatHelp",)

import argparse
import enum
import inspect

from labw_utils.typing_importer import Sequence, Type, Iterable, Callable


class _EnhancedHelpFormatter(argparse.HelpFormatter):
    def _expand_help(self, action: argparse.Action):
        params = {**vars(action), "prog": self._prog}
        for name in list(params):
            if params[name] is argparse.SUPPRESS:
                del params[name]
        for name in list(params):
            if hasattr(params[name], "__name__"):
                params[name] = params[name].__name__  # type: ignore
        if action.choices is not None:
            choices_str = "\nCHOICES:"
            if isinstance(action.choices, enum.EnumMeta):
                for c in action.choices:
                    choices_str += f"\n    {c.name} -- {c.__doc__}"
            elif isinstance(action.choices, Iterable):
                for c in action.choices:
                    choices_str += f"\n    {c}"
            else:
                raise TypeError(f"Type of {action.choices} ({type(action.choices)}) should be Enum or Iterable")
        else:
            choices_str = ""

        help_str = action.help if action.help is not None else ""

        if action.required:
            req_opt_prefix = "[REQUIRED] "
            if action.default is not None and action.default is not argparse.SUPPRESS:
                raise ValueError(f"Argument {action} setted both required and default ({action.default})!")
        else:
            req_opt_prefix = "[OPTIONAL] "
        if not hasattr(action.type, "__name__"):
            dtype_prefix = ""
        elif isinstance(action.type, type):
            dtype_prefix = "Type: " + action.type.__name__ + "; "  # type: ignore
        elif isinstance(action.type, Callable):
            reta = inspect.signature(action.type).return_annotation
            dtype_prefix = "Type: " + str(reta) + "; "  # type: ignore
        else:
            dtype_prefix = "Type: " + action.type.__name__ + "; "  # type: ignore

        default_prefix = ""
        if "%(default)" not in help_str:
            if action.default is not argparse.SUPPRESS:
                if action.required is False:
                    defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                    if action.option_strings or action.nargs in defaulting_nargs:
                        if isinstance(action.default, enum.Enum):
                            default_prefix = f"Default: {action.default.name}"
                        else:
                            default_prefix = f"Default: {action.default}"
                else:
                    default_prefix = "No defaults"
        return (req_opt_prefix + dtype_prefix + default_prefix).strip() + "\n" + help_str % params + choices_str

    def _metavar_formatter(self, action, default_metavar):
        if action.metavar is not None:
            result = action.metavar
        elif action.choices is not None:
            if isinstance(action.choices, enum.EnumMeta):
                choice_strs = [str(choice.name) for choice in action.choices]
            elif isinstance(action.choices, Iterable):
                choice_strs = [str(choice) for choice in action.choices]
            else:
                raise TypeError(f"Type of {action.choices} ({type(action.choices)}) should be Enum or Iterable")

            result = "{%s}" % ",".join(choice_strs)
        else:
            result = default_metavar

        def format(tuple_size):
            if isinstance(result, tuple):
                return result
            else:
                return (result,) * tuple_size

        return format

    def _format_action(self, action: argparse.Action):
        help_position = min(self._action_max_length + 2, self._max_help_position)
        action_header = "%*s" % (self._current_indent, "") + self._format_action_invocation(action) + "\n"
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = help_text.split("\n")
            parts.append("%*s%s\n" % (help_position, "", help_lines[0]))
            for line in help_lines[1:]:
                parts.append("%*s%s\n" % (help_position, "", line))
        elif not action_header.endswith("\n"):
            parts.append("\n")
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


_ACTION_GROUP_TILE_REPLACEMENT_DICT = {
    "optional arguments": "OPTIONS",  # Python < 3.10
    "options": "OPTIONS",  # Python >= 3.10
    "positional arguments": "PARAMETERS",
    None: "",
}


class ArgumentParserWithEnhancedFormatHelp(argparse.ArgumentParser):
    """"""

    def format_help(self) -> str:
        """
        :raises ValueError: If an argument set both ``required`` and ``default``.
        """
        formatter = _EnhancedHelpFormatter(prog=self.prog)
        formatter.add_text(self.description)

        formatter.add_usage(
            usage=self.usage, actions=self._actions, groups=self._mutually_exclusive_groups, prefix="SYNOPSIS: "
        )

        for action_group in self._action_groups:
            formatter.start_section(_ACTION_GROUP_TILE_REPLACEMENT_DICT.get(action_group.title, action_group.title))
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        formatter.add_text(self.epilog)
        return formatter.format_help()

    def to_markdown(self) -> str:
        """Work in process -- do not use."""

        class State:
            current_indent = 0
            current_buff = ""

        s = State()

        def append(line: str):
            if line is None:
                line = ""
            s.current_buff += " " * (s.current_indent * 4)
            s.current_buff += line.strip()
            s.current_buff += "\n"

        def indent():
            s.current_indent += 1

        def dedent():
            if s.current_indent > 0:
                s.current_indent -= 1

        def bold(instr: str) -> str:
            return "**" + instr.replace("*", "\\*") + "**"

        def italic(instr: str) -> str:
            return "*" + instr.replace("*", "\\*") + "*"

        def code(instr: str) -> str:
            assert "`" not in instr
            return "`" + instr + "`"

        append(self.description)
        append("")
        formatter = _EnhancedHelpFormatter(prog=self.prog)
        formatter._width = 10240  # Now you're long enough
        usage = (
            bold("SYNOPSIS:")
            + " "
            + code(
                formatter._format_usage(
                    usage=self.usage,
                    actions=self._actions,
                    groups=self._mutually_exclusive_groups,
                    prefix="",
                ).strip()
            )
        )
        append(usage)
        append("")

        for action_group in self._action_groups:
            if not action_group._group_actions:
                continue
            append(bold(_ACTION_GROUP_TILE_REPLACEMENT_DICT.get(action_group.title, action_group.title) + ":"))
            append("")

            for action in action_group._group_actions:
                ah = action.help if action.help is not None else "NO HELP"

                required_str = italic("[REQUIRED]") if action.required else italic("[OPTIONAL]")
                if not action.option_strings:
                    append("- " + code(action.dest) + " " + required_str + ": " + ah)
                else:
                    append("- " + ", ".join(map(code, action.option_strings)) + " " + required_str + ": " + ah)

                indent()
                append("")
                if not hasattr(action.type, "__name__"):
                    dtype_prefix = ""
                else:
                    dtype_prefix = bold("Type:") + " " + code(action.type.__name__)  # type: ignore
                append(dtype_prefix)

                append("")
                default_prefix = ""
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                    if action.option_strings or action.nargs in defaulting_nargs:
                        default_prefix = bold("Default:") + " " + code(repr(action.default))
                append(default_prefix)

                append("")
                if action.choices is not None:
                    append(bold("Choices:"))
                    append("")
                    for c in action.choices:
                        append("- " + code(repr(c)))

                dedent()

            append("")

        append(self.epilog)
        append("")
        return s.current_buff


def enum_to_choices(src_enum: Type[enum.Enum]) -> Sequence[str]:
    """
    Convert an enum to choices.

    >>> class sample_enum(enum.Enum):
    ...     A = 1
    ...     B = 2
    >>> enum_to_choices(sample_enum)
    ['A', 'B']

    .. versionadded:: 1.0.3
    """
    retl = []
    for v in src_enum:
        retl.append(v.name)
    return retl
