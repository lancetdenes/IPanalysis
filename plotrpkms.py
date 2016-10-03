import numpy as np
import matplotlib.pyplot as plt

#IP always Y input is always X
#
def IP_rpkms_to_logplots(IP_file, input_file):
    IP_num = IP_file[0]
    input_num = input_file[0]
    IP_dict = {}
    input_dict = {}
#parses IP and input rpkm table files into dictionaries      
    for line in open(IP_file):
        vals = line.strip().split('\t')
        genename = vals[0]
        if not genename == '#Gene':
            rpkm = float(vals[2])
            IP_dict[genename] = rpkm

    for line in open(input_file):
        vals2 = line.strip().split('\t')
        genename = vals2[0]
        if genename in IP_dict:
            rpkm = float(vals2[2])
            input_dict[genename] = rpkm

    input = []
    IP = []
    genenames = []
    in_log_list = []
    IP_log_list = []
    toplist_IP = []
    toplist_in = []
    bottomlist_in = []
    bottomlist_IP = []
    toplist = []
    bottomlist = []
    
#writes file with rpkm values, logged rpkm values, relevant ratios, and whether or not they are enriched between IP and input
    file = open('datafile_%s_%s.txt' % (IP_num, input_num), 'w')
    file.write('gene\trpkm_input\trpkm_IP\tin_log\tIP_log\tIPtoin\tIPtoinlog\tside\n')
    for key in input_dict:

        rpkm_input = input_dict[key]
        rpkm_IP = IP_dict[key]        
        input.append(rpkm_input)
        IP.append(rpkm_IP)
        genenames.append(key)
        in_log = np.log10(rpkm_input)
        IP_log = np.log10(rpkm_IP)
        in_log_list.append(in_log)
        IP_log_list.append(IP_log)
        IPtoin = (rpkm_IP/rpkm_input)
        IPtoinlog = (IP_log/in_log)
        #calculate if in top (enriched in IP) or bottom (enriched in input)
        if IP_log > ((0.8*in_log)-0.5):
            toplist_IP.append(IP_log)
            toplist_in.append(in_log)
            toplist.append(key)
            toporbottom = 'top'
        else:
            bottomlist_IP.append(IP_log)
            bottomlist_in.append(in_log)
            bottomlist.append(key)
            toporbottom = 'bottom'
        file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (key, rpkm_input, rpkm_IP, in_log, IP_log, IPtoin, IPtoinlog, toporbottom))
    file.close()
    
#makes plot with top and bottom separated
    plt.figure(1)
    plt.subplot(212)
    plt.plot(bottomlist_in, bottomlist_IP, 'bo')
    plt.subplot(211)
    plt.plot(toplist_in, toplist_IP, 'bo')
    plt.savefig('low_vs_high_%svs%s.jpg' %(IP_num, input_num))
    plt.close()
#writes a file with a list of all genes in the top and all genes in the bottom    
    with open('topgenes_%s_%s.txt' %(IP_num, input_num), 'w') as topgenes:
        for gene in toplist:
            topgenes.write('%s\n' % gene)
   
    
    with open('bottomgenes_%s_%s.txt' %(IP_num, input_num), 'w') as bottomgenes:
        for gene in bottomlist:
            bottomgenes.write('%s\n' % gene)
   


# run on each pair of samples 
'''
IP_rpkms_to_logplots('2.bam.rpkm', '1.bam.rpkm')
IP_rpkms_to_logplots('4.bam.rpkm', '3.bam.rpkm')
IP_rpkms_to_logplots('6.bam.rpkm', '5.bam.rpkm')
IP_rpkms_to_logplots('8.bam.rpkm', '7.bam.rpkm')
'''
