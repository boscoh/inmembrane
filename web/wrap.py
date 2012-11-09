import sys
import os
import markdown


def do(fname):
  text_md = open(fname).read()
  text_html = markdown.markdown(text_md)
  template = open('template.html').read()
  html = template % { 'insert': text_html }
  html_fname = os.path.splitext(fname)[0] + '.html'
  open(html_fname, 'w').write(html)
  os.system('open ' + html_fname)


if __name__ == '__main__':
  [do(f) for f in sys.argv[1:]]
