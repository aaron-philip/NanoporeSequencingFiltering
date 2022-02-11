# Nanopore Sequencing Filtering
Associated files used in work to smooth the raw signals from Nanopore DNA Sequencing with the intent of achieving more accurate base-calling.

The goal of this project is to introduce a flexible structure capable of taking in raw data from Nanopore DNA Sequencing, processing/refining the signals, and 
passing the processed data to the Guppy base-caller. We intend to use filtering techniques inspired by circuit structures. Currently, the data refinement techniques 
present in base-callers like Guppy or other deep-learning identifiers will be well suited to one type of genomic or instrument data but yields poor results when 
applied in other contexts. Our intent is to introduce filtering techniques that can adjust to different types of data to achieve high accuracy reads.

