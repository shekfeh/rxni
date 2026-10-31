"""Structural evidence for audited rearrangement blocks in runtime_v6_7.

Do not import training/audit labels here. Only mapped reactant-product graph
and local structural context may determine event truth. Reagent/conditions are
not consulted. Installation is strictly gated by frozen positive/boundary tests.
"""
from __future__ import annotations
from rdkit import Chem


def _parts(facts):
    rx = facts.get('reaction')
    if not isinstance(rx, str) or rx.count('>') != 2:
        return None
    r, _, p = rx.split('>')
    rm, pm = Chem.MolFromSmiles(r), Chem.MolFromSmiles(p)
    if rm is None or pm is None:
        return None
    return (rm, pm, _mapped(rm), _mapped(pm))


def _mapped(m):
    return {a.GetAtomMapNum(): a for a in m.GetAtoms() if a.GetAtomMapNum() > 0}


def _bond(a, b):
    if a is None or b is None:
        return 0.0
    x = a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
    return float(x.GetBondTypeAsDouble()) if x else 0.0


def _nb(a):
    return {b.GetAtomMapNum() for b in a.GetNeighbors() if b.GetAtomMapNum() > 0}


def _migrated(r, p, roots={'C','Si'}):
    for m in r.keys() & p.keys():
        root = r[m]
        if root.GetSymbol() not in roots or root.GetAtomicNum() != p[m].GetAtomicNum():
            continue
        for src in _nb(root) - _nb(p[m]):
            for tgt in _nb(p[m]) - _nb(root):
                if src in r and src in p and tgt in r and tgt in p:
                    yield m, src, tgt


def _env_conserved(r, p, root, src, tgt):
    # Exclude the two exchanged anchors; other mapped attachments must persist.
    old = sorted((x, _bond(r[root], r[x])) for x in _nb(r[root]) - {src, tgt} if x in r and x in p)
    new = sorted((x, _bond(p[root], p[x])) for x in _nb(p[root]) - {src, tgt} if x in r and x in p)
    return old == new


def _anchors_conserved(r,p,atom,exclude):
    return (_nb(r[atom]) - {exclude}) == (_nb(p[atom]) - {exclude})


def _carbonyl(a):
    return a.GetAtomicNum() == 6 and any(
        n.GetAtomicNum() == 8 and _bond(a, n) >= 1.8 for n in a.GetNeighbors()
    )


def _single_oxygen(a):
    return any(n.GetAtomicNum() == 8 and _bond(a,n) <= 1.2 for n in a.GetNeighbors())


def _element_neighbors(a, elements):
    return any(n.GetSymbol() in elements for n in a.GetNeighbors())


def _rmaps(m, ring):
    return {m.GetAtomWithIdx(i).GetAtomMapNum() for i in ring}


def _event(name, subtype, *maps, center=None, **extra):
    z = {'event': name, 'role': 'specific_event', 'subtype': subtype,
         'support_maps': sorted({int(x) for x in maps if int(x) > 0}),
         'mapping_scope': 'v67_frozen_graph_structural_rearrangement_contract'}
    if center is not None:
        z['center_map'] = int(center)
    z.update(extra)
    return z


def _wolff(rm, pm, r, p):
    out=[]
    for d,rda in r.items():
        if rda.GetAtomicNum() != 6 or d not in p:
            continue
        ns=[n for n in rda.GetNeighbors() if n.GetAtomicNum()==7 and _bond(rda,n)>=1.8
            and any(nn.GetAtomicNum()==7 and _bond(n,nn)>=1.8 for nn in n.GetNeighbors() if nn.GetIdx()!=rda.GetIdx())]
        if not ns or not _carbonyl(p[d]):
            continue
        for ac in rda.GetNeighbors():
            a=ac.GetAtomMapNum()
            if a not in p or not _carbonyl(ac):
                continue
            if not (_bond(rda,ac)>0 and _bond(p[d],p[a])>0):
                continue
            if _carbonyl(p[a]):
                continue
            if any(nn.GetAtomicNum()==7 and _bond(p[d],nn)>=1.8 for nn in p[d].GetNeighbors()):
                continue
            out.append(_event('wolff_rearrangement','alpha_diazoketone_to_ketene_reorganization',a,d,center=d,acyl_c_map=a,diazo_c_map=d))
    return out


def _peroxide_atoms(rm):
    return {a.GetAtomMapNum() for a in rm.GetAtoms()
            if a.GetAtomicNum()==8 and a.GetAtomMapNum()>0
            and any(x.GetAtomicNum()==8 and _bond(a,x)<=1.2 for x in a.GetNeighbors())}


def _criegee_explicit(rm,pm,r,p):
    out=[]
    for mig,s,t in _migrated(r,p,{'C'}):
        if r[s].GetAtomicNum()!=6 or r[t].GetAtomicNum()!=8:
            continue
        if not (_bond(r[s],r[t])>0 and _bond(p[s],p[t])>0):
            continue
        if not any(n.GetAtomicNum()==8 and n.GetIdx()!=r[s].GetIdx()
                   for n in r[t].GetNeighbors()):
            continue
        if not _env_conserved(r,p,mig,s,t):
            continue
        out.append(_event('criegee_rearrangement',
            'explicit_peroxide_c_to_o_migration',mig,s,t,center=s))
    return out


def _baeyer_villiger(rm,pm,r,p):
    out=[]
    peroxide=_peroxide_atoms(rm)
    for a,ac in r.items():
        if a not in p or not (_carbonyl(ac) and _carbonyl(p[a])):
            continue
        # A ketone, not an ester, acid, amide or acyl halide.
        singles=[n for n in ac.GetNeighbors() if _bond(ac,n)<1.8]
        if len([n for n in singles if n.GetAtomicNum()==6])!=2:
            continue
        if any(n.GetAtomicNum() in (7,8,9,17,35,53) for n in singles):
            continue
        for mig in singles:
            m=mig.GetAtomMapNum()
            if m not in p or _bond(p[a],p[m]):
                continue
            for o in p[a].GetNeighbors():
                om=o.GetAtomMapNum()
                if om not in peroxide or o.GetAtomicNum()!=8:
                    continue
                if _bond(p[a],o)>=1.8 or _bond(o,p[m])<0.8:
                    continue
                # Inserted peroxide O must originate from actual R-O-O connectivity.
                if any(_bond(ac,r[om])>0 for _ in [0] if om in r):
                    continue
                out.append(_event('baeyer_villiger_oxygen_insertion',
                    'ketone_to_ester_or_lactone',a,m,om,center=a,
                    migrating_c_map=m,inserted_o_map=om))
    return out


def _ra1_ra2(rm,pm,r,p,facts):
    out=[]
    for m,s,t in _migrated(r,p,{'C'}):
        if r[s].GetAtomicNum()!=6 or r[t].GetAtomicNum()!=6:
            continue
        if not (_bond(r[s],r[t])>0 and _bond(p[s],p[t])>0):
            continue
        if not _env_conserved(r,p,m,s,t):
            continue
        if _single_oxygen(r[s]) and _single_oxygen(r[t]) and not _carbonyl(r[s]) and _carbonyl(p[s]):
            out.append(_event('semipinacol_like_migration','oxygenated_skeletal_c_migration',m,s,t,center=s,migrator_map=m,target_map=t))
            continue
        # Strict RA2: target-side O/X, aromatic migrator, same C framework;
        # exclude new heterocyclic annulation and broad mapped rewirings.
        if (r[m].GetIsAromatic() and r[t].GetAtomicNum()==6 and not r[t].GetIsAromatic()
                and _element_neighbors(r[t],{'O','F','Cl','Br','I'})
                and r[s].GetAtomicNum()==6 and not r[s].GetIsAromatic()
                and not _carbonyl(r[s])
                and not facts.get('ring_formations')
                and len(facts.get('formed_bonds', [])) == 1
                and len(facts.get('broken_bonds', [])) == 1):
            out.append(_event('wagner_meerwein_like_migration','target_activated_c_skeletal_migration',m,s,t,center=t,migrator_map=m,source_map=s))
    return out


def _ra3(rm,pm,r,p):
    out=[]
    for rr in rm.GetRingInfo().AtomRings():
        rset=_rmaps(rm,rr)
        if 0 in rset or len(rr)<3:
            continue
        for pp in pm.GetRingInfo().AtomRings():
            pset=_rmaps(pm,pp)
            if len(pp)!=len(rr)+1 or 0 in pset or not rset.issubset(pset):
                continue
            new=list(pset-rset)
            if len(new)!=1 or new[0] not in r:
                continue
            ex=new[0]
            if r[ex].GetAtomicNum()!=6 or not any(n.GetAtomicNum()==7 for n in r[ex].GetNeighbors()):
                continue
            if any(n.GetAtomicNum()==7 for n in p[ex].GetNeighbors()):
                continue
            for ctr in r[ex].GetNeighbors():
                c=ctr.GetAtomMapNum()
                if c not in rset or c not in p:
                    continue
                if _single_oxygen(r[c]) and _carbonyl(p[c]):
                    out.append(_event('tiffeneau_demjanov_ring_expansion','tiffeneau_like_ring_expansion',ex,c,center=c,reactant_ring_size=len(rr),product_ring_size=len(pp)))
                elif _single_oxygen(p[c]) and not _carbonyl(p[c]):
                    out.append(_event('tiffeneau_demjanov_ring_expansion','demjanov_like_ring_expansion',ex,c,center=c,reactant_ring_size=len(rr),product_ring_size=len(pp)))
    return out


def _ra4(rm,pm,r,p):
    out=[]
    for m,s,t in _migrated(r,p,{'C'}):
        if not (r[m].GetIsAromatic() and _bond(r[s],r[t])>0
                and _carbonyl(r[s]) and _carbonyl(r[t])
                and not _carbonyl(p[t]) and _single_oxygen(p[t])):
            continue
        if not _env_conserved(r,p,m,s,t):
            continue
        # Migrating source becomes carboxylate/carboxylic-acid-like.
        pso=[n for n in p[s].GetNeighbors() if n.GetAtomicNum()==8]
        if _carbonyl(p[s]) and len(pso)>=2:
            out.append(_event('benzilic_acid_like_rearrangement','dicarbonyl_aryl_migration',m,s,t,center=s))
    return out


def _re_rf(rm,pm,r,p):
    out=[]
    for m,s,t in _migrated(r,p,{'C','Si'}):
        if not _env_conserved(r,p,m,s,t):
            continue
        root=r[m];sr=r[s];tr=r[t];sel=sr.GetSymbol();tel=tr.GetSymbol()
        if root.GetSymbol()=='Si' and (sel,tel) in {('C','O'),('O','C')}:
            if _bond(r[s],r[t]) and _bond(p[s],p[t]):
                out.append(_event('brook_retro_brook_rearrangement','silyl_'+sel.lower()+'_to_'+tel.lower(),m,s,t,center=m))
            continue
        if root.GetSymbol()!='C':continue
        if _carbonyl(root) and _carbonyl(p[m]):
            if (sel,tel) in {('O','N'),('N','O')}:
                # Anchors on O and N must survive except acyl-root bond.
                if _anchors_conserved(r,p,s,m) and _anchors_conserved(r,p,t,m):
                    out.append(_event('acyl_heteroatom_relocation','acyl_'+sel.lower()+'_to_'+tel.lower(),m,s,t,center=m))
            if (sel,tel)==('O','C') and not tr.GetIsAromatic():
                # Oxygen stays bound to its non-root group, and target is alpha to a carbonyl.
                if (_anchors_conserved(r,p,s,m) and _anchors_conserved(r,p,t,m)
                        and any(_carbonyl(nn) for nn in tr.GetNeighbors())):
                    aromatic_source=any(nn.GetIsAromatic() and nn.GetAtomicNum()==6 for nn in sr.GetNeighbors())
                    if aromatic_source:
                        subtype='baker_venkataraman_acyl_o_to_c'
                    else:
                        # Enol oxygen must stay bound to an unsaturated carbon.
                        has_enol=any(nn.GetAtomicNum()==6 and any(_bond(nn,x)>=1.8 for x in nn.GetNeighbors() if x.GetAtomicNum()==6) for nn in sr.GetNeighbors())
                        subtype='enol_o_to_c_acyl_shift' if has_enol else None
                    if subtype:
                        out.append(_event('acyl_o_to_c_rearrangement',subtype,m,s,t,center=m))
            if (sel,tel)==('O','C') and tr.GetIsAromatic():
                # Steglich-specific heteroaromatic acyl shift: source acyloxy
                # becomes a carbonyl on its existing ring carbon, and the
                # target aromatic carbon becomes nonaromatic with an acyl bond.
                if (not p[t].GetIsAromatic() and
                    any(nn.GetAtomicNum()==6 and nn.GetIsAromatic() and
                        nn.GetAtomMapNum() in p and _bond(p[s],p[nn.GetAtomMapNum()])>=1.8
                        for nn in sr.GetNeighbors())):
                    out.append(_event('steglich_like_acyl_rearrangement',
                        'heteroaromatic_acyl_o_to_c_shift',m,s,t,center=m))
        if root.GetIsAromatic() and (sel,tel)==('O','N'):
            # Chapman: O and N on the same imidate carbon before migration.
            if (any(nn.GetAtomicNum()==6 and nn.GetIdx()!=root.GetIdx()
                    and _bond(nn,tr)>0 for nn in sr.GetNeighbors())
                and _anchors_conserved(r,p,s,m) and _anchors_conserved(r,p,t,m)):
                out.append(_event('chapman_rearrangement','aryl_o_to_n_imidate_transfer',m,s,t,center=m))
        if (sel,tel) in {('N','C'),('S','C')}:
            if (sr.GetFormalCharge()>0 and tr.GetFormalCharge()<0
                and _bond(sr,tr)>0 and _bond(p[s],p[t])>0):
                out.append(_event('stevens_rearrangement','onium_ylide_alkyl_migration',m,s,t,center=t))
    return out


def _claisen(rm,pm,r,p):
    """Aromatic allyl ether -> ortho-allyl phenol connectivity.

    Deliberately do not demand the pericyclic terminal-C mapping: atom mappers
    may preserve the allyl chain and attach the original O-bound carbon instead.
    The v14 map-free counterfactual independently validates this endpoint class.
    """
    out=[]
    for o,ra in r.items():
        if ra.GetAtomicNum()!=8 or o not in p:continue
        for anchor in [n for n in ra.GetNeighbors() if n.GetAtomicNum()==6 and n.GetIsAromatic()]:
            am=anchor.GetAtomMapNum()
            if am not in p or _bond(ra,anchor)!=1 or _bond(p[o],p[am])<0.8:continue
            if p[o].GetTotalNumHs()==0:continue
            for alpha in ra.GetNeighbors():
                al=alpha.GetAtomMapNum()
                if al not in p or alpha.GetAtomicNum()!=6 or alpha.GetIsAromatic() or _bond(ra,alpha)!=1 or _bond(p[o],p[al]):continue
                for beta in alpha.GetNeighbors():
                    bm=beta.GetAtomMapNum()
                    if bm not in p or beta.GetAtomicNum()!=6 or _bond(alpha,beta)!=1:continue
                    for gamma in beta.GetNeighbors():
                        gm=gamma.GetAtomMapNum()
                        if gm not in p or gamma.GetAtomicNum()!=6 or _bond(beta,gamma)!=2:continue
                        chain={al,bm,gm}
                        if not (_bond(p[al],p[bm])>0 and _bond(p[bm],p[gm])>0):continue
                        if not any(_bond(p[x],p[y])>=1.8 for x,y in [(al,bm),(bm,gm)]):continue
                        for ortho in anchor.GetNeighbors():
                            om=ortho.GetAtomMapNum()
                            if om not in p or ortho.GetAtomicNum()!=6 or not ortho.GetIsAromatic():continue
                            if _bond(anchor,ortho)!=1.5:continue
                            attached=[x for x in chain if _bond(p[om],p[x])>0 and _bond(ortho,r[x])==0]
                            if attached:
                                out.append(_event('claisen_aromatic_rearrangement','ortho_allyl_aryl_ether_to_phenol',o,am,al,bm,gm,om,center=am))
    return out


def _piancatelli(rm,pm,r,p):
    out=[]
    for ring in rm.GetRingInfo().AtomRings():
        if len(ring)!=5:continue
        atoms=[rm.GetAtomWithIdx(i) for i in ring]
        if not (sum(a.GetAtomicNum()==8 for a in atoms)==1 and sum(a.GetAtomicNum()==6 for a in atoms)==4 and all(a.GetIsAromatic() for a in atoms)):
            continue
        furan_m={a.GetAtomMapNum() for a in atoms if a.GetAtomicNum()==6}
        if 0 in furan_m:continue
        for ra in atoms:
            if ra.GetAtomicNum()!=6:continue
            for ctr in ra.GetNeighbors():
                if ctr.GetIdx() in ring or ctr.GetAtomicNum()!=6 or ctr.GetIsAromatic():continue
                cm=ctr.GetAtomMapNum()
                if cm not in p:continue
                if not any(nn.GetAtomicNum()==8 and nn.GetTotalNumHs()>0 for nn in ctr.GetNeighbors()):continue
                aryl=[nn for nn in ctr.GetNeighbors() if nn.GetAtomicNum()==6 and nn.GetIsAromatic() and nn.GetIdx() not in ring]
                if not aryl:continue
                for pp in pm.GetRingInfo().AtomRings():
                    if len(pp)!=5 or _rmaps(pm,pp)!=furan_m|{cm} or not all(pm.GetAtomWithIdx(i).GetAtomicNum()==6 for i in pp):continue
                    ks=[]; hs=[]
                    for idx in pp:
                        a=pm.GetAtomWithIdx(idx)
                        ks.extend([idx for n in a.GetNeighbors() if n.GetIdx() not in pp and n.GetAtomicNum()==8 and _bond(a,n)>1.8])
                        hs.extend([idx for n in a.GetNeighbors() if n.GetIdx() not in pp and n.GetAtomicNum()==8 and _bond(a,n)<1.3 and n.GetTotalNumHs()>0])
                    if not ks or not hs:continue
                    if not any(_bond(p[cm],p[a.GetAtomMapNum()])>0 for a in aryl if a.GetAtomMapNum() in p):continue
                    out.append(_event('piancatelli_like_rearrangement','furylcarbinol_to_hydroxycyclopentenone',cm,*sorted(furan_m),center=cm))
    return out


def detect_v67_rearrangement_extensions(facts):
    """New positive-evidence events; no conditions or reaction-name labels."""
    z=_parts(facts)
    if z is None:return []
    rm,pm,r,p=z
    out=[]
    for fn in (_wolff,_criegee_explicit,_baeyer_villiger,_ra1_ra2,_ra3,_ra4,_re_rf,_claisen,_piancatelli):
        if fn is _ra1_ra2:
            out.extend(fn(rm,pm,r,p,facts))
        else:
            out.extend(fn(rm,pm,r,p))
    # Keep distinct chemistry and mapped supports; collapse duplicates only.
    uniq={}
    for a in out:uniq[(a['event'],a['subtype'],tuple(a['support_maps']))]=a
    return list(uniq.values())

# Explicit primary ownership policy. Valid orthogonal architecture is preserved.
RRG_V67_EVENTS = (
    'wolff_rearrangement', 'criegee_rearrangement',
    'baeyer_villiger_oxygen_insertion', 'semipinacol_like_migration',
    'wagner_meerwein_like_migration', 'tiffeneau_demjanov_ring_expansion',
    'benzilic_acid_like_rearrangement', 'brook_retro_brook_rearrangement',
    'acyl_heteroatom_relocation', 'chapman_rearrangement',
    'acyl_o_to_c_rearrangement', 'steglich_like_acyl_rearrangement',
    'stevens_rearrangement', 'claisen_aromatic_rearrangement',
    'piancatelli_like_rearrangement',
)

RRG_CONFLICTS = {
    'wolff_rearrangement': {'carbonyl_deoxygenation_or_reduction', 'fragmentation', 'multi_bond_fragmentation'},
    'baeyer_villiger_oxygen_insertion': {'oxidative_deformylation', 'o_alkylation', 'carbonyl_substitution'},
    'wagner_meerwein_like_migration': {'elimination'},
    'tiffeneau_demjanov_ring_expansion': {'fragmentation', 'multi_bond_fragmentation'},
    'acyl_heteroatom_relocation': {'n_deprotection','amidation','carbonyl_substitution'},
    'acyl_o_to_c_rearrangement': {'carbonyl_substitution'},
    'steglich_like_acyl_rearrangement': {'carbonyl_substitution'},
    'benzilic_acid_like_rearrangement': {'carbonyl_substitution'},
    'chapman_rearrangement': {'aryl_or_heteroaryl_amination'},
    'stevens_rearrangement': {'n_dealkylation'},
    'piancatelli_like_rearrangement': {'arene_partial_hydrogenation', 'o_dealkylation'},
    'claisen_aromatic_rearrangement': {'o_dealkylation'},
    'carbon_to_nitrogen_migration': {'carbonyl_substitution'},
}


def resolve_v67_rearrangement_primary(chosen):
    """Suppress only competing primary labels supported by an audited owner.

    All original detector observations remain in `specific_events`. Orthogonal
    ring/oxidation/redox events remain eligible. This function is deliberately
    not a generic priority-based eraser.
    """
    owners=[x for x in chosen if x.get('event') in RRG_CONFLICTS]
    if not owners:
        return chosen
    out=[]
    for rec in chosen:
        rname=rec.get('event')
        support=set(rec.get('support_maps') or [])
        center=rec.get('center_map')
        reject=False
        for owner in owners:
            if rname not in RRG_CONFLICTS[owner['event']]:continue
            ours=set(owner.get('support_maps') or [])
            owner_center=owner.get('center_map')
            # Fragmentation uses a different center convention; require
            # a multi-map overlap rather than suppressing a distinct split.
            fragmentation=rname in {'fragmentation','multi_bond_fragmentation'}
            shared=support & ours
            if fragmentation:
                reject=len(shared)>=2
            elif (owner['event']=='wolff_rearrangement'
                  and rname=='carbonyl_deoxygenation_or_reduction'):
                reject=center is not None and int(center) in ours
            else:
                reject=(center is not None and owner_center is not None
                        and int(center)==int(owner_center)) or len(shared)>=2
            if reject:break
        if not reject:out.append(rec)
    return out
