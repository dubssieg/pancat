build:
  @git pull && pip install --upgrade pip && python -m pip install -r requirements.txt --upgrade && python -m pip install .