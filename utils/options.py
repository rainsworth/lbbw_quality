# Options for killms pipeline code

import ConfigParser
import os
import struct
import re

def download_cat(path,url, **kwargs): #need to catch cases where files are packed - what about csv? - conversion?
    if kwargs['cat'].upper() == 'LOTSS':
        import download_lotss_catalogue
        download_lotss_catalogue(kwargs['msin'], ResultsFile=kwargs['lotss_output'], Radius=float(kwargs['lotss_radius']), AllFile=kwargs['lotss_output']+'_all', DoDownload=kwargs['lotss_download'])
    else:
        filename = url.split('/')[-1]
        if filename.split('.')[-1] == 'gz':
            fitsname = '.'.join(filename.split('.')[:-1])
        else:
            fitsname = filename

    path_to_file = path + fitsname
    if os.path.isfile(path_to_file):
        print("Catalogue "+fitsname+" already exists - skipping download")
    else:
        os.system("wget " + url +" "+ path)
        if filename.split('.')[-1] == 'gz':
            os.system("gunzip -c " + filename + " > "+ path_to_file)
    return path + fitsname

def _get_terminal_size_linux():
    ''' From https://gist.github.com/jtriley/1108174 '''
    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import termios
            cr = struct.unpack('hh',
                               fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
            return cr
        except:
            pass
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            return None
    return int(cr[1]), int(cr[0])

def options(optlist,option_list):

    # option_list format is: section, name, type, default
    # section names are used in the output dict only if names are not unique

    odict = {}
    config=ConfigParser.SafeConfigParser()
    filenames=[]
    cmdlineset=[]
    if isinstance(optlist,str):
        optlist=[optlist]

    for o in optlist:
        if o[:2]=='--':
            optstring=o[2:]
            result=re.match('(\w*)-(\w*)\s*=\s*(.*)',optstring)
            if result is None:
                print 'Cannot parse option',optstring
            else:
                cmdlineset.append(result.groups())
        else:
            filenames.append(o)

    config.read(filenames)
    for c in cmdlineset:
        try:
            config.add_section(c[0])
        except ConfigParser.DuplicateSectionError:
            pass
        config.set(c[0],c[1],c[2])
    cased={int: config.getint, float: config.getfloat, bool: config.getboolean, str: config.get, list: lambda x,y: eval(config.get(x,y))}
    for o in option_list:
        (section, name, otype, default)=o[:4]
        # if this name is duplicated in another section, we need to know
        count=0
        for o2 in option_list:
            if o2[1]==name: count+=1
        # get result
        try:
            result=cased[otype](section,name)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            result=default
        if count>1:
            odict[section+'_'+name]=result
        else:
            odict[name]=result
    if odict['logging']=='None':
        odict['logging']=None
    return odict

def typename(s):
    return str(s).replace("type '","").replace("'","")

def print_options(option_list):
    import textwrap
    from auxcodes import bcolors
    # expected to be called if a config file is not specified. Print a
    # list of options
    width,height=_get_terminal_size_linux()
    sections=sorted(set(x[0] for x in option_list))
    klen=max([len(x[1]) for x in option_list])
    tlen=max([len(typename(x[2])) for x in option_list])
    fstring='%-'+str(klen)+'s = %-'+str(tlen)+'s (default %s)'
    indent=' '*(klen+3)
    for s in sections:
        print bcolors.OKBLUE+'\n[%s]' % s+bcolors.ENDC
        for o in option_list:
            if len(o)==4:
                (section, name, otype, default)=o
                doc=None
            elif len(o)==5:
                (section, name, otype, default, doc)=o
            else:
                print 'Oops!',o
                continue
            if section==s:

                print bcolors.BOLD+fstring % (name, typename(otype), str(default))+bcolors.ENDC
                if doc is not None:
                    print textwrap.fill(doc,width-1,initial_indent=indent,subsequent_indent=indent)

if __name__=='__main__':
    import sys
    from parset import option_list
    config=sys.argv[1:]
    print options(config,option_list)
