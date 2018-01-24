# New release packaging checklist:

* Always use a clean checkout, otherwise setup.py may pick up files not
  tracked by git that you don't actually want to distribute.

* Make sure CHANGELOG is updated with major changes since the last 
  release (look through the commit history)

* Update version in inmembrane/__init__.py.

* git commit -a

* Create a new tag for the release:
  git tag inmembrane-0.xx
  git push --tags
  git push

* Run:
  * virtualenv /tmp/inmembrane_venv
  * source /tmp/inmembrane_venv/bin/activate
  * python setup.py sdist
  * pip uninstall inmembrane
    pip install dist/inmembrane-<version>.tar.gz
  * Test the installed version: 
    inmembrane_scan --test

As per https://packaging.python.org/tutorials/distributing-packages/#upload-your-distributions

* Create a `~/.pypirc` file like:

```
[pypi]
username=<my_username>
password=<my_password>
~                           
```

* chmod 600 ~/.pypirc

* Push a new version to PyPi:
  * pip install twine
  * twine upload dist/*

* Switch to the the gh-pages branch:
  git checkout gh-pages

* Update the pypi download link in index.html.wrap. 
  Regenerate HTML docs with:
  python wrap.py *html.wrap
  
* Commit and push the changes.
  git commit -a; git push
