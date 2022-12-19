FROM bentoml/model-server:0.11.0-py37
LABEL author="ersilia"

RUN pip install crem@git+https://github.com/DrrDom/crem
RUN pip install 'rdkit>=2017.09'

WORKDIR /repo
COPY . /repo
