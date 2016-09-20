import glob
import setuptools

setuptools.setup(name='frea',
                 version='0.7.2',
                 description='Functional Region Enrichment Analysis',
                 long_description='Tools for interpreting GWAS using regulatory annotations',
                 author='Abhishek Sarkar',
                 author_email='aksarkar@mit.edu',
                 url='http://web.mit.edu/aksarkar/',
                 packages=setuptools.find_packages(),
                 scripts=glob.glob('bin/*'),
                 entry_points={
                     'console_scripts': [
                         'ucsc_to_impg=frea.summary.process:ucsc_to_impg',
                         'diagram=frea.summary.process:diagram',
                         'cardiogram=frea.summary.process:cardiogram',
                         'iibdgc=frea.summary.process:iibdgc',
                         'frea-build-table=frea.matched:build_table',
                         'frea-output-pheno=frea.simulation:output_oxstats_pheno',
                         'frea-output-summary=frea.simulation:output_summary',
                         'frea-impg-maps=frea.impg:oxstats_legend_to_impg_map',
                         'frea-convert-haps=frea.impg:oxstats_haps_to_impg_haps',
                         'frea-split-haps=frea.impg:split_impg_haps',
                         'frea-impg-holdout=frea.summary.process:scz_to_impg'
                     ]
                 }
                 )
