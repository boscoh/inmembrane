import sys
import os
import markdown
import re


usage = """
Converts HTML template files (using .wrap) extension
with a simple tag {{tag.md}} for markdown file
includes.

Usage: python wrap.py *html.wrap
"""


def do(wrap_fname):
  html = open(wrap_fname).read()
  for match in re.findall(r'{{.*}}', html):
    fname = match[2:-2]
    if os.path.isfile(fname):
      text_md = open(fname).read()
    else:
      text_md = ""
    text_html = '\n' + markdown.markdown(text_md) + '\n'
    html = html.replace(match, text_html)
  html_fname = wrap_fname.replace('.wrap', '')
  open(html_fname, 'w').write(html)
  print "Made", os.path.abspath(html_fname)


if __name__ == '__main__':
  if len(sys.argv) < 2:
    print usage
  else:
    for f in sys.argv[1:]:
      if f.endswith('wrap'):
        do(f)
      else:
        print "Processes html templates with .wrap extension"
