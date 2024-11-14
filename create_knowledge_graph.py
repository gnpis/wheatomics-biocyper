from biocypher import BioCypher, Resource
from template_package.adapters.gene_adapter import (
    WheatomicsAdapter,
    WheatomicsAdapterNodeType,
    WheatomicsAdapterEdgeType,
    WheatomicsAdapterGeneField,
    WheatomicsAdapterTranscriptField,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher()

# Download and cache resources (change the directory in the options if needed)
# urls = "https://file-examples.com/wp-content/storage/2017/02/file_example_CSV_5000.csv"
# resource = Resource(
#     name="Example resource",  # Name of the resource
#     url_s=urls,  # URL to the resource(s)
#     lifetime=7,  # seven days cache lifetime
# )
# paths = bc.download(resource)  # Downloads to '.cache' by default
# print(paths)
# You can use the list of paths returned to read the resource into your adapter

# # Choose node types to include in the knowledge graph.
# # These are defined in the adapter (`adapter.py`).
# node_types = [
#     WheatomicsAdapterNodeType.GENE,
#     NodeType.ACCESSION,
#     NodeType.TRANSCRIPT,
#     NodeType.CONDITION,
#     NodeType.ACCESSION,
# ]

# # Choose protein adapter fields to include in the knowledge graph.
# # These are defined in the adapter (`adapter.py`).
# node_fields = [
#     # Gene
#     GeneField.NAME,
#     GeneField.LOCUS,
#     GeneField.ACCESSION_ID,
#     GeneField.GENOME_VERSION,
#     GeneField.ANONOTATION_VERSION,
#     #Accession
#     AccessionField.NAME,
#     # Transcript
#     TranscriptField.POSITION,
#     TranscriptField.ABUNDANCE_VALUE,
#     TranscriptField.GENE_TRANSCRIPT_RELATION, # REALLY??
#     TranscriptField.CONDITION,
    


# # Transcripts
# # - position
# # - abundance value
# # - gene-transcript-rel
# # - condition
# # - metadata.tsv
# # 	- development_stage
# # 	- Organ
# # 	- Traitement (stress)


#     # Diseases
#     ExampleAdapterDiseaseField.ID,
#     ExampleAdapterDiseaseField.NAME,
#     ExampleAdapterDiseaseField.DESCRIPTION,
# ]

# # Choose protein adapter fields to include in the knowledge graph.
# # These are defined in the adapter (`adapter.py`).
# node_fields = [
#     # Proteins
#     ExampleAdapterProteinField.ID,
#     ExampleAdapterProteinField.SEQUENCE,
#     ExampleAdapterProteinField.DESCRIPTION,
#     ExampleAdapterProteinField.TAXON,
#     # Diseases
#     ExampleAdapterDiseaseField.ID,
#     ExampleAdapterDiseaseField.NAME,
#     ExampleAdapterDiseaseField.DESCRIPTION,
# ]

# edge_types = [
#     ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
#     ExampleAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION,
# ]

# Create a protein adapter instance
# adapter = WheatomicsAdapterNodeType(
adapter = WheatomicsAdapter(
    # node_types=node_types,
    # node_fields=node_fields,
    # edge_types=edge_types,
    # we can leave edge fields empty, defaulting to all fields in the adapter
)


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

bc.write_schema_info(as_node=True)

# Print summary
# bc.summary()
bc.log_missing_input_labels()
