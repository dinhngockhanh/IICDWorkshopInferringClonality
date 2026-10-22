vb"""Point Python 2.7's Modules/Setup SSL block at /opt/openssl11."""

from pathlib import Path

setup = Path("Modules/Setup.dist")
text = setup.read_text()
old = (
    "#SSL=/usr/local/ssl\n"
    "#_ssl _ssl.c \\\n"
    "#\t-DUSE_SSL -I$(SSL)/include -I$(SSL)/include/openssl \\\n"
    "#\t-L$(SSL)/lib -lssl -lcrypto\n"
)
new = (
    "SSL=/opt/openssl11\n"
    "_ssl _ssl.c \\\n"
    "\t-DUSE_SSL -I$(SSL)/include -I$(SSL)/include/openssl \\\n"
    "\t-L$(SSL)/lib -Wl,-rpath,/opt/openssl11/lib -lssl -lcrypto\n"
)
if old not in text:
    raise SystemExit("Could not find the commented SSL block in Modules/Setup.dist")
setup.write_text(text.replace(old, new, 1))
Path("Modules/Setup").write_text(setup.read_text())
