"""Read the documented cascade tables (.docx) and return event-level chains.

Only the standard library is used. A .docx is a zip; the first table in
word/document.xml is read, the ID column is forward-filled (as in the paper's
R code) and the 'Proposed Encoding' column is parsed into ordered peril pairs.
"""
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import OrderedDict

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def read_first_table(path):
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    tbl = root.find(".//" + W + "tbl")
    rows = []
    for tr in tbl.findall(W + "tr"):
        cells = []
        for tc in tr.findall(W + "tc"):
            cells.append("".join(t.text or "" for t in tc.iter(W + "t")).strip())
        rows.append(cells)
    return rows


def parse_encoding(enc):
    """'SEQ,TS;TS,CFL;' -> [('SEQ','TS'), ('TS','CFL')] (order kept, duplicates kept)."""
    pairs = []
    for seg in enc.split(";"):
        pr = [p.strip() for p in seg.split(",")]
        if len(pr) >= 2 and pr[0] and pr[1]:
            pairs.append((pr[0], pr[1]))
    return pairs


def load_chains(path):
    """Return OrderedDict event_id -> list of chains (each chain = list of pairs)."""
    rows = read_first_table(path)
    header = [h.lower() for h in rows[0]]
    enc_col = next(i for i, h in enumerate(header) if "encoding" in h)
    id_col = next(i for i, h in enumerate(header) if h.startswith("id") or h == "id")
    events = OrderedDict()
    cur = None
    for r in rows[1:]:
        if len(r) <= enc_col:
            continue
        raw_id = r[id_col].strip() if id_col < len(r) else ""
        # In the mainland tables the ID cell of one event is sometimes split over
        # two rows, the second holding only a region label (东/南/西/北/中) or an
        # 'and <id>' continuation. Those rows belong to the preceding event.
        is_continuation = bool(re.fullmatch(r"[东南西北中]", raw_id)) or raw_id.startswith("and ")
        if raw_id and not is_continuation:
            cur = raw_id
        enc = r[enc_col].strip()
        if not enc:
            continue
        events.setdefault(cur, []).append(parse_encoding(enc))
    return events


def load_many(paths):
    out = OrderedDict()
    for p in paths:
        for k, v in load_chains(p).items():
            out.setdefault(k, []).extend(v)
    return out


def event_edge_sets(events, drop_self_loops=True):
    """event_id -> set of directed edges documented in that event."""
    es = OrderedDict()
    for ev, chains in events.items():
        s = set()
        for ch in chains:
            for a, b in ch:
                if drop_self_loops and a == b:
                    continue
                s.add((a, b))
        es[ev] = s
    return es


def event_node_sets(events):
    ns = OrderedDict()
    for ev, chains in events.items():
        s = set()
        for ch in chains:
            for a, b in ch:
                s.add(a); s.add(b)
        ns[ev] = s
    return ns
