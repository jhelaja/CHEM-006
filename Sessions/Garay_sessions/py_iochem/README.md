# py-ioChem
Diego Garay-Ruiz, 2024 -
This library contains utilities to interface ioChem-BD's REST API and automate the
retrieval of calculations from the platform and the input/output of reaction networks

## Automated retrieval of calculations
Calculations in ioChem-BD are organized in *collections*, which can be either public (Find and Browse modules) or private (Create module).

*py-ioChem* enables the user to:
1. Explore the structure of the collection -> calculation names and contained files.
2. Automatically retrieve all *input* files and CML *outputs*.
3. Filter the query to include or exclude calculations containing given substrings, to avoid downloading too many files for large collections.

### Find and Browse modules
For these *public* collections, two main pieces of input are required:
- The *handle* of the calculation, present in the URL of the collection -> a code of the form:
    - **Find**: 10/XXXXXX
    - **Browse**: 100/XXXXX
- The URL for the REST API:
    - **Find**: `https://www.iochem-bd.org/rest`.
    - **Browse**: depends on the instance, which might or might have not it available. 
        For the BSC instance, `https://iochem-bd.bsc.es/rest` -> contact your administrator for other instances.

With this, the basic procedure to get *all* files in the collection is:

```python 
from py_iochem import RESTAPIManager
rest_url = "https://www.iochem-bd.org/rest"
handle = "10/XXXXXX"                                       # Find-like, but it is analogous for Browse

rest_manager = RESTAPIManager.CollectionHandler(rest_url,handle,token=None,service="Find")
rest_manager.get_items()                                   # Necessary to get the structure of the collection
track_info,_ = rest_manager.get_files(outdir="cmlfiles")   # Directory where input and output files will be saved
rest_manager.save_tracking(track_info,output_file="tracking_list.csv")
```

From this code snippet, all files will be saved to the *cmlfiles* directory, while a *tracking_list.csv* file will be
generated containing tabulated data for each file: calculation name, calculation ID, output file name and original file name in 
the ioChem-BD collection. An example would be:
```
calculation1,AAAA,input1.in,input1.in
calculation1,AAAA,cml_AAAA.cml,output.cml
calculation1_test,BBBB,input1_test.in,input1_test.in
calculation1_test,BBBB,cml_BBBB.cml,output.cml
calculation2,CCCC,input2.in,input2.in
calculation2,CCCC,cml_CCCC.cml,output.cml
calculation2_freq,DDDD,input2_freq.in,input2_freq.in
calculation2_freq,DDDD,cml_DDDD.cml,output.cml
```
where the *input* files preserve their original names and the output files are named depending on the handle of the calculation.

To **filter** the batch download, the user should prepare a list of the substrings that should be filtered. 
For instance: 
```python 
from py_iochem import RESTAPIManager
rest_url = "https://www.iochem-bd.org/rest"
handle = "10/XXXXXX"                                       # Find-like, but it is analogous for Browse

filter_words = ["test","freq"]                             # Look for these words in the calculation names
rest_manager = RESTAPIManager.CollectionHandler(rest_url,handle,token=None,service="Find")
rest_manager.get_items()                                   # Necessary to get the structure of the collection
rest_manager.filter_names(filter_words,"exclude")          # Calculations matching the name filter will NOT be downloaded
track_info,_ = rest_manager.get_files(outdir="cmlfiles")   # Directory where input and output files will be saved
rest_manager.save_tracking(track_info,output_file="tracking_list.csv")
```

And then, all files with "test" or "freq" in their names will be excluding, only downloading:
```
calculation1,AAAA,input1.in,input1.in
calculation1,AAAA,cml_AAAA.cml,output.cml
calculation2,CCCC,input2.in,input2.in
calculation2,CCCC,cml_CCCC.cml,output.cml
```

If we applied the same snippet but with:
```python
rest_manager.filter_names(filter_words,"include")          # ONLY calculations matching the name filter WILL BE downloaded
```
we would get:
```
calculation1_test,BBBB,input1_test.in,input1_test.in
calculation1_test,BBBB,cml_BBBB.cml,output.cml
calculation2_freq,DDDD,input2_freq.in,input2_freq.in
calculation2_freq,DDDD,cml_DDDD.cml,output.cml
```

It is advisable to explore the collection first, either using the ioChem-BD GUI or via Python, before proceeding to download
a large amount of files, to properly use the files.

To explore the calculations in a Python script, the user should check the elements contained in the `rest_manager.itemList` attribute:
a list of dictionaries containing properties for all the calculations.

### Create module
For these *private* collections, while the procedure is mostly analogous to the one outlined in the previous section,
there are four parts of input.
 - The *id* of the calculation, which can be checked in the ioChem-BD main interface: collection properties, lower right part of the screen, as **Id XXXXX**
 - The URL for the REST API, which depends on the instance, which might or might have not it available. 
   For the BSC instance, `https://iochem-bd.bsc.es/create/rest` -> contact your administrator for other instances.
 - The *username* (e-mail address) of the current user, used for logging in.
 - The *password* for the user -> this will be passed privately, via the *getpass* module, and not shown as plain text anywhere.

 As **Create** is a private module, credentials are necessary to retrieve data. Thus, the generation of the API manager here requires identification, as follows:
 ```python
    from py_iochem import RESTAPIManager,prompt_token
    create_url = "https://iochem-bd.bsc.es/create/rest/"
    create_user = "foo@bar.es"
    create_token = prompt_token(create_url,create_user)
    rest_manager = RESTAPIManager.CollectionHandler(create_url,create_id,create_token,"Create")
 ``` 
 where the *prompt_token* function will ask the user for interactive password entry: if credentials are right, the bearer token required for the API identification will be generated,
 and the *rest_manager* object will work just like before -> the same filtering functions are available.

## Reaction network format
ioChem-BD manages graphs in JSON format, compliant with Cytoscape. 

To map nodes and edges in the graphs to calculations, there must be a `cid_formula` field where for every node/edge,
the identifiers in ioChem-BD (`calcId`) of the calculations used to define that node/edge are specified and separated by *+* symbols.

For example, for a reaction A + B --TS1--> C:

ioChem-BD collection:
- A		cid = 12
- B		cid = 13
- C		cid = 14
- TS1	cid = 15

*Node 1*. `formula`: A+B, `cid_formula`: 12+13  
*Node 2*. `formula`: C, `cid_formula`: 14  
*Edge 1 - 2*. `formula`: TS1, `cid_formula`: 15

Therefore, to upload an arbitrary reaction network to ioChem-BD, three steps must be taken:
1. Generation of JSON-formatted graph, with `formula` fields for all nodes/edges.
2. Automated upload of all calculations to ioChem-BD
3. Generation of Reaction Graph report 

*py-ioChem* contains utilities to manage this process, using only a configuration file to 
define input parameters (user name, report name, graph file...). The upload part requires a 
shell script relying on ioChem-BD's shell client, that should be installed separately.

