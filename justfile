build:
  @python setup.py install

dist:
  @git add * && git commit -m "Just auto-commit" && git push