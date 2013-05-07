from subprocess import Popen, PIPE;
class Viterbi:
	def __init__(self, extraArgs = ""):
		self.process = Popen(
			['perl', 'viterbi.pl'],
			cwd='.',
			shell=0,
			stdin=PIPE, stdout=PIPE)

	def tag(self, input):
		self.process.stdin.write(input + '\n')
		self.process.stdin.flush()
		return self.process.stdout.readline()

v = Viterbi()
