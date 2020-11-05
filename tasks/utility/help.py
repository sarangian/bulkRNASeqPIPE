import luigi


def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output

class manual(luigi.Task):
	#Global Parameters
	def requires(self):
		return []

	def output(self):
		return []

	def run(self):
		module="ENA"
		print(run_cmd(module))