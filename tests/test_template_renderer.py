from autogenu import template_renderer


def test_generated_file_uses_utf8_and_lf_line_endings(monkeypatch, tmp_path):
    template_directory = tmp_path / "templates"
    template_directory.mkdir()
    (template_directory / "sample.in").write_bytes(b"value={{value}}\r\n")
    monkeypatch.setattr(template_renderer, "TEMPLATE_DIRECTORY", template_directory)
    output_path = tmp_path / "generated.txt"

    template_renderer.write_generated_file(output_path, "sample.in", value=42)

    assert output_path.read_bytes() == b"value=42\n"
