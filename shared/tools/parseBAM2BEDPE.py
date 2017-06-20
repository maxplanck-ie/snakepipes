import enum
import os.path
import sys

class bam2BEDPE:
    STRING_PIPE = '|'
    STRING_REDIRECT = '>'
    BEDPE_SUFFIX = 'bedpe'
    BAM_SUFFIX = 'bam'

    def __init__(self, sourcefile, do_adjust = True, samtools_command = 'samtools'):
        self.flags = ""
        self.doInclude = True
        self.do_adjust = do_adjust

        self.paths = enum.Enum
        self.paths.SAMTOOLS = samtools_command

        self.sourcefile = sourcefile
        self.targetfile = '.'.join([self.sourcefile.replace('.'+BAM_SUFFIX,''),BEDPE_SUFFIX])

        self.offsetUpstream = 0
        self.offsetDownstream = 0

    def _prepare_decipher_bam(self):
        if not self.fileExists():
            raise IOError("\n".join(["File not found:", self.sourcefile]))

        return(' '.join([self.paths.SAMTOOLS, 'view', self.flags, self.sourcefile]))

    def _prepare_extract(self):
        awk_chrom = '$3'
        awk_start = '$4 - 1 + ' + str(self.offsetUpstream)
        awk_end = '$4 + abs($9)  - 1 + ' + str(self.offsetDownstream)

        awk_func = "function abs(v) {return v < 0 ? -v : v}"
        awk_print = "".join(['\'', awk_func, '{','print','(',', '.join([awk_chrom, awk_start, awk_end]),')' '}','\''])

        return(' '.join(["awk","-v OFS='\\t'", awk_print]))

    def _prepare_target(self):
        return(str(self.targetfile))

    def set_path_to_samtools(self, path):
        self.paths.SAMTOOLS = path

    def fileExists(self):
        return(os.path.isfile(self.sourcefile))

    def _setATAC_offset(self):
        print("Setting ATAC cutsite offset")
        self.offsetUpstream = -4
        self.offsetDownstream = +5

    def selectFirstMateOnly(self):
        self.flags = '-f 0x40'

    def do_adjust_cutsites(self):
        self.do_adjust = True
        self._setATAC_offset()

    def compileCommand(self):
        self.command = ' '.join([   self._prepare_decipher_bam(), bam2BEDPE.STRING_PIPE, \
                                    self._prepare_extract(), bam2BEDPE.STRING_REDIRECT, \
                                                                self._prepare_target()])
    def execute(self):
        import subprocess
        print('Executing:\n',self.command)

        if len(self.command) > 0:
            process = subprocess.Popen(self.command, stdout=subprocess.PIPE, shell = True)
            self.output, self.error = process.communicate()
        else:
            IOError("Nothing to execute", self.command)

    def __str__(self):
        return("".join(['Prepared for execution','\n', '\t', self.command]))
