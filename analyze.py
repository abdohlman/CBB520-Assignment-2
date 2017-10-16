#!/usr/bin/python
import os, sys, shutil
import cgitb, cgi
import time, datetime
from subprocess import Popen, PIPE
#import subprocess as sub

cgitb.enable()

# Initial parameters
BARCODE = None
DEFAULT = 'SRR4841864'
TRUNCATE_TF = True #False
KMERSIZE = 50

# Define directories
STKDIR = '/home/abd30/sratoolkit/bin/'
TMPDIR = '/usr/lib/cgi-bin/tmp/'
ABYDIR = '/usr/lib/cgi-bin/assembly/'

# Locate reference genome
REFGENOME = '/usr/lib/cgi-bin/genome/S288C.fna'

def tag(text,tag='p'):
        '''Helper function for returning text flanked by html tag. Default to <p>.'''
        return '<{0}>{1}</{0}>'.format(tag, text)

def display(msg='',tag='p',tstep=0.5):
        '''Helper function for displaying messages to client in real time.'''
        print '<{0}>{1}</{0}>'.format(tag, msg)
        sys.stdout.flush()
        time.sleep(tstep)

def execute(cmd):
        '''Helper function to quickly execute shell commands and return stdout and stderr.'''
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        output, error = p.communicate()
        return output, error

def wait_for(p, tstep=10):
        '''Wait for a process <p> to finish, printing every <tstep> seconds if still running'''
        t1 = time.time()
        sys.stdout.write('\r>> Please wait.')
        sys.stdout.flush()
        while p.poll() is None:
                time.sleep(tstep)
                sys.stdout.write('\r.')
                sys.stdout.flush()
                t2 = time.time()
        hours, rem = divmod(t2-t1,3600)
        minutes, seconds = divmod(rem,60)
        display('>> Took {0} minutes and {1} seconds to complete.'.format(int(minutes),int(seconds)))
        return

def get_file(truncate=TRUNCATE_TF):
        '''Get accession barcode from index.html and download it. If none is received, use the default. Return fastq path.'''
        form = cgi.FieldStorage()
        usrinput = form.getvalue('query') 
	truncate = form.getvalue('truncate')
	global BARCODE
	try:
                BARCODE = form.getvalue('query')
		if BARCODE[:3] != 'SRR':
			raise ValueError()
                display(tag("Getting {0}...".format(BARCODE),'b'))
        except:
                BARCODE = DEFAULT
                display(tag("Unrecognized barcode provided, using {0} as default.".format(DEFAULT),'b'))
        if truncate == 'on':
		display('>> Truncating input to 100k lines.')
                cmd = [STKDIR+'fastq-dump','--split-spot','-X','100000',BARCODE, '-O', TMPDIR ]
        else:
                cmd = [STKDIR+'fastq-dump','--split-spot',BARCODE,'-O',TMPDIR ]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, bufsize=0)
        wait_for(p)
        output, error = p.communicate()
	line1, line2 = output.split('Written') #.replace(BARCODE,BARCODE+'.')
        display('>> '+line1[:-1]+'.')
	display('>> Wrote'+line2[:-1]+'.')
        display('>> Short read archive SRR4841864 was downloaded successfully.')
        fastqpath = TMPDIR + BARCODE + '.fastq'
        return fastqpath

def pair_good(score1,score2):
        '''Return True if pair has >50 bases of quality score >25 criteria, else False'''
        values1, values2 = zip(*[ (ord(x1)-33,ord(x2)-33) for x1,x2 in zip(score1,score2) ])
        values1, values2 = list(values1), list(values2)
        if min(values1 + values2) >= 25:
                return True
        for lst in [values1, values2]:
                bad_ix = [0] + [ ix for ix in range(len(lst)) if lst[ix] < 25 ] + [101]
                gapsizes = [  abs(bad_ix[i-1]-bad_ix[i]) for i in range(1,len(bad_ix)) ]
                if max(gapsizes) < 50:
                        return False
        return True

def count_quality_pairs(fastqpath):
        '''Read a fastq path line by line and determine if both pairs have >50 bases of quality score >25.'''
        display("<b>Counting quality pairs...</b>")
        output, error = execute(['wc','-l',fastqpath])
        nlines = int(output.split(' ',1)[0])
        if nlines % 8 != 0:
                display('>> WARNING: {0} may be corrupted! nlines % 8 = {1}'.format(fastqpath,str(nlines % 8)))
        nreads = nlines/8
        ct = 0
        t0 = time.time()
        with open(fastqpath,'r') as f:
                t1 = time.time()
                sys.stdout.write('\r>> Please wait.')
                sys.stdout.flush()
                while True:
                        t2 = time.time()
                        if t2 - t1 > 10.0:
                                print '\r.'
                                sys.stdout.flush()
                                t1 = time.time()
                        lines = [ f.readline().rstrip() for x in range(8) ]
                        if not lines[-1]: break
                        score1, score2 = lines[3], lines[7]
                        if pair_good(score1, score2): ct += 1
        hours, rem = divmod(t2-t0,3600)
        minutes, seconds = divmod(rem,60)
        display('>> Took {0} minutes and {1} seconds to complete.'.format(int(minutes),int(seconds)))
        pct = round(100*ct/float(nreads),2)
        display(">> Found a total of {0} quality pairs out of {1} ({2}%)".format(ct, nreads, pct))
        return

def align_to_reference(fastqpath):
        '''Execute BWA mem and pipe to samtools view to output bam file containing alignment output. Return bam path.'''
        display("<b>Aligning to reference genome S288c...</b>")
	bampath = TMPDIR + BARCODE + '.bam'
        with open(bampath,'w') as out:
                p1 = Popen(['bwa','mem','-p',REFGENOME,fastqpath], stdout=PIPE, stderr=PIPE)
                p2 = Popen(['samtools','view','-S','-b','-'], stdin=p1.stdout, stdout=out, stderr=PIPE)
                wait_for(p2)
                for line in p1.stderr.read().split('\n')[:-1]:
                        display('>> '+line)
                display('>> '+p2.stderr.read())
        display('>> Converted fastq to bam.')
        return bampath

def sort_bam(bampath):
        '''Sort a bam file. Return sorted bam path.'''
        display('>> Sorting bam file.')
        sbampath = TMPDIR + BARCODE + '.sorted.bam'
        cmd = ['samtools','sort',bampath,'-f',sbampath]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        wait_for(p)
        return sbampath

def call_snps(sbampath):
        '''Use samtools mpileup and bcftools view to call SNPS. Return vcf path'''
        display("<b>Calling SNPs...</b>")
        vcfpath = TMPDIR + BARCODE + '.vcf'
        with open(vcfpath,'w') as out:
                p1 = Popen(['samtools','mpileup','-B','-u','-f',REFGENOME,'-o',TMPDIR,sbampath], stdout=PIPE, stderr=PIPE)
                p2 = Popen(['bcftools','view','-vcgI','-'], stdin=p1.stdout, stdout=out, stderr=PIPE)
                wait_for(p2)
                display('>> '+p1.stderr.read())
                display('>> '+p2.stderr.read())
        return vcfpath

def count_snps(vcfpath):
        '''Parse vcf file and count number of each type of SNPs. eg. A -> T. Return 4x4 table.'''
        bases = ['A','T','G','C']
        table = [ [0] * 4 for i in range(4) ]
        with open(vcfpath,'r') as f:
                t1 = time.time()
                sys.stdout.write('\r>> Please wait.')
                sys.stdout.flush()
                while True:
                        t2 = time.time()
                        if t2 - t1 > 10:
                                sys.stdout.write('\r.')
                        sys.stdout.flush()
                        t1 = time.time()
                        line = f.readline().rstrip()
                        if not line: break
                        elif line[0] == '#': continue
                        else: ref, alt = line.split('\t')[3:5]
                        for snp in alt.split(','):
                                i,j = bases.index(ref), bases.index(snp)
                                table[i][j] += 1
        hours, rem = divmod(t2-t1,3600)
        minutes, seconds = divmod(rem,60)
        display('>> Took {0} minutes and {1} seconds to complete.'.format(int(minutes),int(seconds)))
        nsnps = sum(sum(row) for row in table)
        display('>> Identified a total of '+str(nsnps)+' SNPs.')
        return table

def show_snps(table):
        '''Take a 4x4 table of SNP counts and nicely display it in html.'''
        display('>> SNP Breakdown (rows: reference base, columns: sequence base)')
        print '<table cellpadding="10">'
        print '<tr>'
        bases = ['A','T','G','C']
        for base in [''] + bases: print tag(base,'th')
        print '</tr>'
        for i in range(4):
                print '<tr>'
                print tag(bases[i],'th')
                for item in table[i]:
                        print tag(str(item),'td')
                print '</tr>'
        print '</table>'

def identify_unmapped(sbampath):
        '''Read sorted bam file and use samtools view to flag unmapped reads. Return unmapped bam file.'''
        display(tag('Calling unmapped reads...','b'))
        ubampath = TMPDIR + BARCODE + '.unmapped.bam'
        with open(ubampath,'w') as out:
                p = Popen(['samtools','view','-b','-f','4',sbampath], stdout=out, stderr=PIPE)
                wait_for(p)
                display('>> [samtools] Saved unmapped reads to '+ubampath)
        return ubampath

def assemble_denovo(ubampath, k=64):
        '''Read unmapped bam file and performs ABYSS de novo assembly on unmapped reads. Return unmapped fasta path.'''
        display(tag('Performing de novo assembly on unmapped reads...','b'))
        display('>> Writing to '+ABYDIR)
        fastapath = ABYDIR + BARCODE + '.fa'
        cmd = ['abyss-pe','-C',ABYDIR,'name='+BARCODE,'k='+str(k),"in='{0}'".format(ubampath)]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        wait_for(p)
        output, error = p.communicate()
        for line in error.splitlines():
                display('>> '+line)
        display('>> Wrote unmapped reads to {0}'.format(fastapath))
        return fastapath

if __name__ == "__main__":

        # Dont litter
        for f in os.listdir('/usr/lib/cgi-bin/assembly/'):
                output, error = execute(['rm',ABYDIR + f])
        for f in os.listdir('/usr/lib/cgi-bin/tmp/'):
                output, error = execute(['rm',TMPDIR + f])

        # Begin html document
        sys.stdout.write('Content-Type: text/html;charset=utf-8\r\n\r\n')
        print
        print '<html>'
        print '<header align="center" style="font-family: helvetica">'
        print '<h1>Assignment 2 Pipeline</h1>'
        print tag('Created by Anders Dohlman for CBB520','p')
        print tag('Duke University 2017','p')
        print '</header>'
        print '<body style="font-family: courier">'

        # Execute pipeline
        fastqpath = get_file()
        count_quality_pairs(fastqpath)
        bampath = align_to_reference(fastqpath)
        sbampath = sort_bam(bampath)
        vcfpath = call_snps(sbampath)
        table = count_snps(vcfpath)
        show_snps(table)
        ubampath = identify_unmapped(sbampath)
        fastapath = assemble_denovo(ubampath,k=KMERSIZE)
        display(tag('Pipeline has completed!','b'))

        # End html document
        print """\
        </body></html>
        """
        
                                                                  
                                                                                                        
