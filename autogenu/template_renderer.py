import re
from pathlib import Path

TEMPLATE_DIRECTORY = Path(__file__).with_name("templates")
PLACEHOLDER = re.compile(r"{{([A-Za-z_][A-Za-z0-9_]*)}}")


def render_template(template_name, **values):
    """Render a packaged text template using explicit ``{{name}}`` placeholders."""
    template = (TEMPLATE_DIRECTORY / template_name).read_text(encoding="utf-8")
    required = set(PLACEHOLDER.findall(template))
    missing = required.difference(values)
    if missing:
        names = ", ".join(sorted(missing))
        raise KeyError(f"Missing template values for {template_name}: {names}")

    rendered = PLACEHOLDER.sub(lambda match: str(values[match.group(1)]), template)
    unresolved = PLACEHOLDER.findall(rendered)
    if unresolved:
        names = ", ".join(sorted(set(unresolved)))
        raise ValueError(f"Unresolved placeholders in {template_name}: {names}")
    return rendered


def write_generated_file(path, template_name, **values):
    """Render a template with stable UTF-8 and LF output."""
    output = render_template(template_name, **values)
    with Path(path).open("w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(output)
