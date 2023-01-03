# CReM: Chemically Reasonable Mutations Framework for Structure Generation

## Model identifiers

- Slug: crem-structure-generation
- Ersilia ID: eos4q1a
- Tags: generative, fragment-based, rule-based

## Model description

The framework is an open source implementation of fragment based generative approaches for exploring chemical space while ensuring chemical validity. This framework utilizes a database of known compounds to come up with interchangable fragments based on the context radius of an input molecule to generate new molecules. The generated molecules vary across inputs by order of 1000s. For sake of simplicity, for such molecules generating over 100 molecules, only 100 diverse molecules are returned in the output. This selection is done using mini batch K Means clustering with 100 clusters and molecules closest to the centroid in each cluster are returned.  

- Input: SMILES
- Output: SMILES
- Model type: Rule-based generative model

## Source code

This framework has been published by Polishchuk, P. CReM: chemically reasonable mutations framework for structure generation. J Cheminform 12, 28 (2020). DOI: [https://doi.org/10.1186/s13321-020-00431-w](https://doi.org/10.1186/s13321-020-00431-w)

- [Code](https://github.com/DrrDom/crem) and [Documentation](https://crem.readthedocs.io/en/latest/)
- Fragment database: replacements02_sc2.db - [database](lhttp://www.qsar4u.com/pages/crem.php) created from ChEMBL v22 structures which contain only organic atoms (C,N,O,S,P,F,Cl,Br,I) and have maximum synthetic complexity score (SCScore) 2.

## License

The BSD 3-Clause License applies to all parts of this repository.

## History

Model was incorporateed on December 19, 2022.

## About us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission or [volunteer](https://www.ersilia.io/volunteer) with us!
