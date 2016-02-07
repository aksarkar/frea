import glob
import setuptools

setuptools.setup(name='frea',
                 version='0.2',
                 description='Functional Region Enrichment Analysis',
                 long_description='Tools for interpreting GWAS using regulatory annotations',
                 author='Abhishek Sarkar',
                 author_email='aksarkar@mit.edu',
                 url='http://web.mit.edu/aksarkar/',
                 packages=setuptools.find_packages(),
                 scripts=glob.glob('bin/*'),
                 entry_points={
                     'console_scripts': ['ucsc_to_impg=frea.summary.process:ucsc_to_impg',
                                         'diagram=frea.summary.process:diagram',
                                         'cardiogram=frea.summary.process:cardiogram',
                                         'iibdgc=frea.summary.process:iibdgc']
                 }
                 )
