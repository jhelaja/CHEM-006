from setuptools import setup,find_packages

setup(name='py_iochem',
      version='0.33',
      author="Diego Garay-Ruiz",
      author_email="dgaray@iciq.es",
      description="Python functions to manage ioChem-BD results, including CML file processing, DOT-formatted graph I/O and REST API management",
      packages=["py_iochem"],
      install_requires=['lxml','networkx','requests','saxonche','importlib_resources'],
      include_package_data=True,
      package_data = {"py_iochem":["stylesheets/*.xsl"]}
)
 
