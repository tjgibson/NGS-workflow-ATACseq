import subprocess as sp
import sys
from itertools import product
from snakemake.shell import shell

db = snakemake.params.database.lower()
species = snakemake.params.species.lower()
release = int(snakemake.params.release)
build = snakemake.params.build


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

spec = ("{build}" if int(release) > 75 else "{build}.{release}").format(
    build=build, release=release
)

if db == "ensembl":
	branch = ""
	if release >= 81 and build == "GRCh37":
		# use the special grch37 branch for new releases
		branch = "grch37/"

elif db == "ensembl_bacteria":
    url_base = "ftp://ftp.ensemblgenomes.org/pub/release-{release}/bacteria/"
elif db == "ensembl_fungi":
    url_base = "ftp://ftp.ensemblgenomes.org/pub/release-{release}/fungi/"
elif db == "ensembl_metazoa":
    url_base = "ftp://ftp.ensemblgenomes.org/pub/release-{release}/metazoa/"
elif db == "ensembl_plants":
    url_base = "ftp://ftp.ensemblgenomes.org/pub/release-{release}/plants/"
elif db == "ensembl_protists":
    url_base = "ftp://ftp.ensemblgenomes.org/pub/release-{release}/protists/"

else:
    raise ValueError("invalid database, must be one of ensembl, ensembl_bacteria, ensembl_fungi, ensembl_metazoa, ensembl_plants, ensembl_protists")

suffixes = ""
datatype = snakemake.params.get("datatype", "")
chromosome = snakemake.params.get("chromosome", "")
if datatype == "dna":
    if chromosome:
        suffixes = ["dna.chromosome.{}.fa.gz".format(chromosome)]
    else:
        suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
elif datatype == "cdna":
    suffixes = ["cdna.all.fa.gz"]
elif datatype == "cds":
    suffixes = ["cds.all.fa.gz"]
elif datatype == "ncrna":
    suffixes = ["ncrna.fa.gz"]
elif datatype == "pep":
    suffixes = ["pep.all.fa.gz"]
else:
    raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

if chromosome:
    if not datatype == "dna":
        raise ValueError(
            "invalid datatype, to select a single chromosome the datatype must be dna"
        )

success = False

for suffix in suffixes:
    if db == "ensembl":
		url = "ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{suffix}".format(
			release=release,
			species=species,
			datatype=datatype,
			spec=spec.format(build=build, release=release),
			suffix=suffix,
			species_cap=species.capitalize(),
			branch=branch,
		)
	else:
		url = "{url_base}/fasta/{species}/{datatype}/{species_cap}.{spec}.{suffix}".format(
			url_base=url_base,
			release=release,
			species=species,
			datatype=datatype,
			spec=spec.format(build=build, release=release),
			suffix=suffix,
			species_cap=species.capitalize(),
			branch=branch,
		)

    try:
        shell("curl -sSf {url} > /dev/null 2> /dev/null")
    except sp.CalledProcessError:
        continue

    shell("(curl -L {url} | gzip -d > {snakemake.output[0]}) {log}")
    success = True
    break

if not success:
    print(
        "Unable to download requested sequence data from Ensembl. "
        "Did you check that this combination of species, build, and release is actually provided?",
        file=sys.stderr,
    )
    exit(1)