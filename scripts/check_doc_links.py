#!/usr/bin/env python3
"""Check repository-local Markdown links and embedded images."""
import re
from pathlib import Path
from urllib.parse import unquote

root = Path(__file__).resolve().parent.parent
files = [*root.glob('*.md'), *root.glob('docs/**/*.md'),
         *root.glob('developer/**/*.md'), *root.glob('examples/**/*.md')]
errors = []
count = 0
for source in files:
    text = source.read_text(encoding='utf-8')
    urls = re.findall(r'\]\(([^)]+)\)', text)
    urls += re.findall(r'<img\s+[^>]*?src="([^"]+)"', text)
    for url in urls:
        url = url.split()[0].strip('<>')
        if re.match(r'^(?:https?://|mailto:)', url):
            continue
        path, _, fragment = url.partition('#')
        target = source.parent / unquote(path) if path else source
        count += 1
        if not target.exists():
            errors.append(f'{source.relative_to(root)}: missing {url}')
        elif fragment and target.suffix.lower() == '.md':
            heads = re.findall(r'^#{1,6}\s+(.+?)\s*#*$', target.read_text(encoding='utf-8'), re.M)
            slugs = [re.sub(r'[^\w\- ]', '', re.sub(r'<[^>]*>|[`*]', '', h).lower()).replace(' ', '-') for h in heads]
            if fragment not in slugs:
                errors.append(f'{source.relative_to(root)}: missing anchor {url}')
for error in errors:
    print(error)
print(f'Checked {count} local links; {len(errors)} errors')
raise SystemExit(bool(errors))
