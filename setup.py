#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from setuptools import setup
from stat import S_IEXEC
from glob import glob
from os.path import join
from os import (chmod, rename, stat, devnull, getcwd, environ, chdir)
import sys
from urllib import FancyURLopener
from subprocess import call, Popen, PIPE
from tempfile import mkdtemp
from shutil import rmtree, copy


__version__ = "0.0.1-dev"


# heavily based on lib.util.download_file from github.com/qiime/qiime-deploy
class URLOpener(FancyURLopener):
    def http_error_default(self, url, fp, errcode, errmsg, headers):
        raise IOError(
            'Could not download %s\nPlease ensure the URL is valid and that '
            'you have an active Internet connection.' % url)


def status(msg):
    """ Write message immediately to stdout

        Parameters
        ----------
        msg : string
            message to print to stdout
    """
    sys.stdout.write(msg)
    sys.stdout.write('\n')
    sys.stdout.flush()


# heavily based on lib.util.download_file from github.com/qiime/qiime-deploy
def download_file(URL, dest_dir, local_file, num_retries=4):
    """ General file downloader

        Parameters
        ----------
        URL : string
            url to download the file from
        dest_dir : string
            directory where you want to download the file
        local_file : string
            output filename of the download
        num_retries: integer
            number of times the function will try to
            download the file

        Returns
        -------
        return_code : integer
            exit status for the download 0 = success,
            1 = fail
    """
    status('  Downloading %s...' % local_file)

    url_opener = URLOpener()
    localFP = join(dest_dir, local_file)
    tmpDownloadFP = '%s.part' % localFP

    return_code = 1
    while num_retries > 0:
        try:
            tmpLocalFP, headers = url_opener.retrieve(URL, tmpDownloadFP)
            rename(tmpDownloadFP, localFP)
            return_code = 0
        except IOError:
            if num_retries == 1:
                status('  Download of %s failed.' % URL)
            else:
                status('  Download failed. Trying again... %d tries remain.' %
                       (num_retries - 1))
            num_retries -= 1
        else:
            num_retries = 0
            status('  %s downloaded successfully.' % local_file)
    return return_code


def app_available(app_name):
    """Check if a binary is available and on the user Path

    Parameters
    ----------
    app_name : string
        Name of the binary, i. e. 'ls', 'gcc' etc.

    Returns
    -------
    output : boolean
        False if the binary is not found, True if the binary is found
    """
    # redirect all output to /dev/null so nothing is seen on screen
    devnull_fd = open(devnull, 'w')
    output = True

    try:
        call([app_name], stdout=devnull_fd, stderr=devnull_fd)
    except OSError:
        output = False
    finally:
        devnull_fd.close()

    return output


def system_call(cmd, error_msg):
    """Call `cmd` and return whether it was successful or not.
       This function is taken and modified from qcli (previously
       `qcli_system_call`).

    Parameters
    ----------
    cmd : string
        command for system call
    error_msg : string
        error message output on an unsuccessful
        system call

    """
    proc = Popen(cmd,
                 shell=True,
                 universal_newlines=True,
                 stdout=PIPE,
                 stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    success = return_value == 0
    if not success:
        status("Unable to %s:" % error_msg)
        status("  stdout:")
        for line in stdout.split('\n'):
            status("    " + line)
        status("  stderr:")
        for line in stderr.split('\n'):
            status("    " + line)
    return success


def download_VSEARCH():
    """ Download the VSEARCH executable and set it
        to the scripts directory
    """
    status("Installing VSEARCH...")

    if sys.platform in ['darwin', 'macos']:
        URL = ('https://github.com/torognes/vsearch/releases/download/'
               'v1.1.1/vsearch-1.1.1-osx-x86_64')
    elif sys.platform in ['linux', 'linux2']:
        URL = ('https://github.com/torognes/vsearch/releases/download/'
               'v1.1.1/vsearch-1.1.1-linux-x86_64')
    else:
        status("Platform %r not supported by VSEARCH.\n" % sys.platform)
        return

    return_value = download_file(URL, 'scripts/', 'vsearch')

    # make the file an executable file
    if not return_value:
        chmod('scripts/vsearch', stat('scripts/vsearch').st_mode | S_IEXEC)
        status("VSEARCH installed.\n")
    else:
        status("VSEARCH could not be installed.\n")


def build_SortMeRNA():
    """Download and build SortMeRNA then copy it to the scripts directory"""
    status("Building SortMeRNA...")

    # SortMeRNA's configure script doesn't correctly guess the C/C++ compilers
    # to use on OS X. Try to figure that out here.
    if sys.platform.lower() in ['darwin', 'macos']:
        cxx = 'clang++'
        cc = 'clang'
    elif sys.platform.lower() in ['linux', 'linux2']:
        cxx = 'g++'
        cc = 'gcc'
    else:
        status("Unknown or unsupported platform %r, so cannot "
               "build SortMeRNA.\n" % sys.platform)
        return

    for compiler in cxx, cc:
        if not app_available(compiler):
            status("%r not installed, so cannot build SortMeRNA.\n" % compiler)
            return

    cwd = getcwd()
    scripts = join(cwd, 'scripts')
    cxx_old = environ.get('CXX', '')
    cc_old = environ.get('CC', '')

    try:
        tempdir = mkdtemp()
        if download_file(
            'ftp://ftp.microbio.me/pub/sortmerna-2.0-no-db.tar.gz', tempdir,
                'sortmerna-2.0-no-db.tar.gz'):
                    status("Could not download SortMeRNA, so cannot "
                           "install it.\n")
                    return

        chdir(tempdir)

        if not system_call('tar xzf sortmerna-2.0-no-db.tar.gz',
                           'extract SortMeRNA archive'):
            return

        chdir('sortmerna-2.0')

        environ['CXX'] = cxx
        environ['CC'] = cc

        if not system_call('bash build.sh', 'build SortMeRNA'):
            return

        copy('sortmerna', scripts)
        copy('indexdb_rna', scripts)
        status("SortMeRNA built.\n")
    finally:
        environ['CXX'] = cxx_old
        environ['CC'] = cc_old

        # remove the source
        rmtree(tempdir)
        chdir(cwd)


def catch_install_errors(install_function, name):
    try:
        install_function()
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        exception_type, exception_value = sys.exc_info()[:2]
        status(
            "Skipping installation of %s due to failure while downloading, "
            "building, or installing:\n  %s: %s\n" %
            (name, exception_type.__name__, exception_value))

# don't build any of the non-Python dependencies if the following modes are
# invoked
if all([e not in sys.argv for e in 'egg_info', 'sdist', 'register']):
    catch_install_errors(download_VSEARCH, 'VSEARCH')
    catch_install_errors(build_SortMeRNA, 'SortMeRNA')

classes = """
    Development Status :: 2 - Pre-Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

long_description = ("Deblur: a greedy deconvolution algorithm based on known "
                    "read error profiles")

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='deblur',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Deblur',
      author="Deblur development team",
      author_email="amnonim@gmail.com",
      url='https://github.com/biocore/deblur',
      test_suite='nose.collector',
      packages=['deblur'],
      package_data={},
      scripts=glob('scripts/*'),
      extras_require={'test': ["nose >= 0.10.1", "pep8"],
                      'doc': ["Sphinx >= 1.2.2", "sphinx-bootstrap-theme"]},
      install_requires=['click', 'numpy >= 1.7',
                        'scikit-bio >= 0.2.2, < 0.3.0',
                        'biom-format >= 2.1.3, < 2.2.0',
                        'burrito < 1.0.0',
                        'burrito-fillings == 0.1.0-dev'],
      dependency_links=[('https://github.com/biocore/burrito-fillings/archive/'
                        'master.zip#egg=burrito-fillings-0.1.0-dev')],
      classifiers=classifiers
      )
