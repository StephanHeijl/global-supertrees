# This bash file downloads Uniprot -> TaxID mappings from Uniprot.
# The resulting file is about 3.7GB in plain text

#curl ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz | gunzip | grep NCBI_TaxID
curl ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.2015_03.gz | gunzip | grep NCBI_TaxID