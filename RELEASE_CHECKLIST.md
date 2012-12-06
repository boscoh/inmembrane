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
  * python setup.py sdist
  * sudo pip uninstall inmembrane
    sudo pip install dist/inmembrane-<version>.tar.gz
  * Test the installed version: inmembrane_scan --test
  * python setup.py register

* Push a new version to PyPi:
  * python setup.py sdist upload

* Switch to the the gh-pages branch (git checkout gh-pages), update the pypi download link in index.md, commit and push the changes.

* Switch to the the gh-pages branch (git checkout gh-pages), update the pypi download link in index.md, commit and push the changes.

* Commit and push.
