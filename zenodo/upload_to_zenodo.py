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

r = upload_large_file("iga-meta.tsv.gz")
r = upload_large_file("igg-meta.tsv.gz")
r = upload_large_file("igm-meta.tsv.gz")
r = upload_large_file("epic-iga.tsv.gz")
r = upload_large_file("epic-igg.tsv.gz")
r = upload_large_file("epic-igm.tsv.gz")
r = upload_large_file("README.md")
