import sys
import os
import markdown
import re

def do(template_fname):
  html = open(template_fname).read()
  for match in re.findall(r'{{.*}}', html):
    fname = match[2:-2]
    if os.path.isfile(fname):
      text_md = open(fname).read()
    else:
      text_md = ""
    text_html = '\n' + markdown.markdown(text_md) + '\n'
    html = html.replace(match, text_html)
  html_fname = template_fname.replace('.template', '')
  open(html_fname, 'w').write(html)


if __name__ == '__main__':
  [do(f) for f in sys.argv[1:]]
