for w in dir(workflow):
   print(w)
   print(getattr(workflow, w))
   print("#"*80)

##print(os.path.dirname(workflow.snakefile))
##print(workflow.overwrite_configfile)
