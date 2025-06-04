from py_iochem import RXioChem as rxioc


config_file = "Report_Config.ini"

rxioc.auto_report_from_json("iochem_graph.json","calcid_tracking.dat",config_file)
