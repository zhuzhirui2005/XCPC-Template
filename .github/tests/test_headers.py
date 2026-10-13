import contextlib
import importlib.util
import io
from pathlib import Path
import shutil
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('headers', Path(__file__).parents[1] / 'scripts/headers.py')
headers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(headers)


class HeadersTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding='utf-8')
        return path

    def test_utf8_bom_and_newlines(self):
        path = self.root / '中文.hpp'
        path.write_bytes(b'\xef\xbb\xbf// test\r\nint n;\r')
        self.assertEqual(headers.read_header(path), '// test\nint n;\n')

    def test_gb18030(self):
        path = self.root / '中文.hpp'
        path.write_bytes('// 中文\r\n'.encode('gb18030'))
        self.assertEqual(headers.read_header(path), '// 中文\n')

    def test_invalid_encoding_fails(self):
        path = self.root / 'bad.hpp'
        path.write_bytes(b'\xff')
        with self.assertRaises(UnicodeDecodeError):
            headers.read_header(path)

    def test_no_headers_fails(self):
        with self.assertRaises(ValueError):
            headers.markdown(self.root)

    def test_markdown_fence_tabs_and_stability(self):
        path = self.write('中文.hpp', '// ```\n\tint n;\n')
        original = path.read_bytes()
        result = headers.markdown(self.root)
        self.assertIn('````cpp\n', result)
        self.assertIn('    int n;', result)
        self.assertEqual(result, headers.markdown(self.root))
        self.assertEqual(path.read_bytes(), original)

    def test_pdf_id_changes_with_content(self):
        path = self.write('sample.hpp', 'int n;')
        first = headers.markdown(self.root)
        path.write_text('int m;', encoding='utf-8')
        second = headers.markdown(self.root)
        first_id = next(line for line in first.splitlines() if 'pdf:trailerid' in line)
        second_id = next(line for line in second.splitlines() if 'pdf:trailerid' in line)
        self.assertNotEqual(first_id, second_id)

    def test_missing_header_fails(self):
        self.write('sample.hpp', 'int n;')
        with self.assertRaisesRegex(ValueError, 'Header.hpp'):
            headers.prepare(self.root, self.root)

    def test_missing_dependency_fails(self):
        self.write('Header.hpp', '')
        self.write('LCT.hpp', '')
        target = self.root / 'prepared'
        target.mkdir()
        with self.assertRaisesRegex(ValueError, 'modint.hpp'):
            headers.prepare(self.root, target)

    @unittest.skipUnless(shutil.which('g++'), 'requires g++')
    def test_compilation_preserves_sources_and_dependencies(self):
        self.write('Header.hpp', 'using ll = long long;\n')
        self.write('modint.hpp', 'struct Mint { int x; };\n')
        self.write('LCT.hpp', 'Mint value;\n')
        original = {p.name: p.read_bytes() for p in self.root.glob('*.hpp')}
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(headers.compile_headers(self.root, 'g++', 2), 0)
        self.assertEqual(original, {p.name: p.read_bytes() for p in self.root.glob('*.hpp')})
        self.assertFalse(list(self.root.glob('*.cpp')))

    @unittest.skipUnless(shutil.which('g++'), 'requires g++')
    def test_failure_does_not_hide_other_headers(self):
        self.write('Header.hpp', '')
        self.write('broken.hpp', 'invalid c++ syntax;')
        self.write('valid.hpp', 'int n;')
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(headers.compile_headers(self.root, 'g++', 2), 1)
        self.assertIn('broken.hpp: FAIL', output.getvalue())
        self.assertIn('valid.hpp: PASS', output.getvalue())


if __name__ == '__main__':
    unittest.main()
