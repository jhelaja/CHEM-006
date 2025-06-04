import py_iochem.RXioChem as rxioc

graph_file = "iochem_graph.json"
config_file = "Report_Config.ini"

rxioc.prepare_network_config(graph_file,config_file)
