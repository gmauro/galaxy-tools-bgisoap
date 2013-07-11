galaxy-tools-bgisoap
====================

Wrapping of the BGI SOAP package for use in Galaxy.

--

Feb 5, 2013
Pending bugs:
1. soap1
Full parameter setting takes forever to finish
2. soapsnp
Full parameter setting not working
3. soapdenovo2
Full parameter setting not working
4. soapdenovo-trans
"reads only used for scaffold assembly" not working
5. scaff submodule
Full parameter setting not working
6. soapcoverage
not working
7. SOAPSNV package
Could not review the output files.

--

Backing up data for migration to a new server

No tool is available to move data, workflows, etc to a new server. This would
be useful for when we want to wipe the Postgresql database but keep all of the
published workflows, published pages and shared libraries.

It looks like the metadata is stored in the postgresql database and actual files
are kepted in the database/files folder.

Let's create a minimal Galaxy instance containing data, workflows, pages which
we can use for backups, refreshing a Galaxy instance, testing, etc.
