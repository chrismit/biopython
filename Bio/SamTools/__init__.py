from multiprocessing import Pool
import subprocess, itertools

def _run_cmd(args):
    child = subprocess.Popen(args, stdout=subprocess.PIPE)
    child.wait()
    out = child.communicate()[0]
    return out

class SamTools:
    """A threaded daemon for querying bamfiles
    
    This is a simple interface to allow for easy multithreaded
    querying of a bam file using samtools
    """
    def __init__(self, bam_source, threads=1, binary='samtools'):
        from multiprocessing import Pool
        self.pool = Pool(processes=threads)
        self.binary = binary
        self.threads = threads
        if isinstance(bam_source,list):
            self.bam_source = ' '.join(bam_source)
        else:
            self.bam_source = bam_source
    
    def _group_iter(self, args, chunks):
        #http://stackoverflow.com/questions/7306522/combining-itertools-and-multiprocessing
        it = iter(args)
        def take():
            while 1: yield list(itertools.islice(it,chunks))
        return iter(take().next,[])
        
    def mpileup(self, **kwrds):
        """A call to mpileup in samtools
        
        This will query a bamfile with the available parameters
        offered by samtools mpileup.  Available options to
        mpileup should be passed as keywords.
        Optional keywords:
        callback: a function to send results to
        Special keyword enhancements:
        The 'r' or 'l' option can be provided a list to query
        many positions/bed files at once.
        
        Example:
        from Bio.SamTools import SamTools
        sTools = '/home/chris/bin/samtools'
        hg19 = '/media/chris/ChrisSSD/ref/human/hg19.fa'
        bamSource = '/media/chris/ChrisSSD/accepted_hits.bam'
        st = SamTools(bamSource,binary=sTools,threads=30)
        print st.mpileup(f=hg19,l='/home/chris/tbed.bed') #tbed is a test bed file
        
        #now with a callback, which is advisable to use to process data as it is generated
        def processPileup(pileup):
            print 'to process',pileup
        
        st.mpileup(f=hg19,r=['chr1:%d-%d'%(i,i+1) for i in xrange(2000001,2001001)],callback=processPileup)
        print st.mpileup(f=hg19,r=['chr1:%d-%d'%(i,i+1) for i in xrange(2000001,2000101)])
        """
        special_keywords = set(['r','l','callback'])
        args = [self.binary,'mpileup']
        to_loop = False
        for i in kwrds:
            if i in special_keywords and isinstance(kwrds[i],list):
                to_loop = kwrds[i]
                to_loop_key = ['-%s'%i]
            elif i != 'callback':
                args+=['-%s'%i,kwrds[i]]
        if 'callback' in kwrds:
            cback = kwrds['callback']
        else:
            cback = None
        results = []
        if to_loop:
            for i in to_loop:
                to_send = [args+to_loop_key+[i]+[self.bam_source]]
                if not cback:
                    results.append(self.pool.apply_async(_run_cmd, to_send))
                else:
                    self.pool.apply_async(_run_cmd, to_send,callback=cback)
        else:
            args.append(self.bam_source)
            if not cback:
                results.append(self.pool.apply_async(_run_cmd, [args]))
            else:
                self.pool.apply_async(_run_cmd, to_send,callback=cback)
        self.pool.close()
        self.pool.join()
        if not cback:
            return [i.get() for i in results]
        else:
            return None
#from Bio.SamTools import SamTools
#sTools = '/home/chris/bin/samtools'
#hg19 = '/media/chris/ChrisSSD/ref/human/hg19.fa'
#bamSource = '/media/chris/ChrisSSD/TH1Alignment/NK/accepted_hits.bam'
#st = SamTools(bamSource,binary=sTools,threads=30)
##print st.mpileup(f=hg19,l='/home/chris/tbed.bed') #tbed is a test bed file
#
##now with a callback, which is advisable to use to process data as it is generated
#def processPileup(pileup):
#    print 'to process',pileup
#
#st.mpileup(f=hg19,r=['chr1:%d-%d'%(i,i+1) for i in xrange(2000001,2001001)],callback=processPileup)
##print st.mpileup(f=hg19,r=['chr1:%d-%d'%(i,i+1) for i in xrange(2000001,2000101)])