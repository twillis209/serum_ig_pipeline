import os
import requests

TOKEN = os.environ.get("ZENODO_TOKEN")
deposition_id = 17010403

headers = {
    "Authorization": f"Bearer {TOKEN}",
}

r = requests.get(
    f"https://zenodo.org/api/deposit/depositions/{deposition_id}", headers=headers
)

bucket_url = r.json()["links"]["bucket"]

def upload_large_file(file_path):
    with open(file_path, "rb") as fp:
        r = requests.put(
            "%s/%s" % (bucket_url, os.path.basename(file_path)),
            data=fp,
            headers=headers,
        )
    return r

files_to_upload = [
    "asthma.tsv.gz",
    "pbc.tsv.gz",
    "psc.tsv.gz",
    "ra.tsv.gz",
    "sle.tsv.gz",
    "crohns.tsv.gz",
    "t1d.tsv.gz",
    "uc.tsv.gz",
    "ms.tsv.gz",
    "derm.tsv.gz",
    "hypothy.tsv.gz",
    "celiac.tsv.gz",
    "igan.tsv.gz",
    "lymphocyte-counts.tsv.gz",
    "dennis-iga.tsv.gz",
    "dennis-igg.tsv.gz",
    "eldjarn-iga.tsv.gz",
    "eldjarn-igg.tsv.gz",
    "eldjarn-igm.tsv.gz",
    "liu-iga.tsv.gz",
    "pietzner-iga.tsv.gz",
    "pietzner-igg.tsv.gz",
    "pietzner-igm.tsv.gz",
    "scepanovic-iga.tsv.gz",
    "scepanovic-igg.tsv.gz",
    "scepanovic-igm.tsv.gz",
    "iga-meta.tsv.gz",
    "igg-meta.tsv.gz",
    "igm-meta.tsv.gz",
    "epic-iga.tsv.gz",
    "epic-igg.tsv.gz",
    "epic-igm.tsv.gz",
    "README.md"
]

for file in files_to_upload:
    r = upload_large_file(file)
    print(f"Uploaded {file}: {r.status_code}")
