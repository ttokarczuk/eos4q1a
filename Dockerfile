FROM bentoml/model-server:0.11.0-py37
LABEL author="ersilia"

RUN pip install crem@git+https://github.com/DrrDom/crem.git@6ddfba45ab904f1a4b53610b7d38dd96fd82a1e5
RUN pip install 'rdkit>=2017.09'
RUN pip install scikit-learn==1.0.2

WORKDIR /repo
COPY . /repo
