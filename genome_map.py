#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
import urllib.parse

def parse_args():
    parser = argparse.ArgumentParser(description="Plota genoma com domínios Pfam (Alta Qualidade).")
    parser.add_argument("-g", "--gff", required=True, help="Arquivo GFF3.")
    parser.add_argument("-p", "--pfam", required=True, help="Tabela Pfam/InterProScan.")
    parser.add_argument("-o", "--output", default="genome_plot.pdf", help="Saída (.pdf ou .svg).")
    parser.add_argument("-t", "--title", default="Organização Genômica", help="Título.")
    parser.add_argument("--dpi", type=int, default=300, help="DPI.")
    return parser.parse_args()

def parse_gff_robust(gff_file):
    genes = []
    genome_length = 0
    print(f"[1/4] Lendo GFF (Modo Robusto): {gff_file}...")
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 8: continue # Pula linhas quebradas
            
            feature_type = parts[2]
            if feature_type not in ["gene", "CDS", "mRNA"]: continue
                
            try:
                start, end = int(parts[3]), int(parts[4])
                strand = 1 if parts[6] == "+" else (-1 if parts[6] == "-" else 0)
                
                # Parse de atributos
                attrs = {}
                for attr in parts[8].split(";"):
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attrs[k.strip()] = urllib.parse.unquote(v.strip())
                
                # Tenta pegar o ID mais confiável
                # Prioridade: Name > ID > locus_tag > gene_id gerado
                gene_id = attrs.get("ID", f"gene_{line_num}")
                name = attrs.get("Name", gene_id) 
                
                genes.append({
                    "id": gene_id,
                    "name": name,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
                if end > genome_length: genome_length = end
            except ValueError: continue

    return pd.DataFrame(genes), genome_length

def main():
    args = parse_args()
    
    # 1. Carregar GFF
    genes_df, genome_length = parse_gff_robust(args.gff)
    if genes_df.empty:
        print("ERRO: Nenhum gene encontrado no GFF.")
        return
    print(f"      -> {len(genes_df)} genes encontrados.")

    # 2. Carregar Pfam (Tratando arquivo sem cabeçalho)
    print(f"[2/4] Lendo Pfam: {args.pfam}...")
    try:
        # Define nomes das colunas padrão do InterProScan
        col_names = [
            'seq_id', 'md5', 'len', 'analysis', 'signature_acc', 'domain_name',
            'start_aa', 'end_aa', 'score', 'status', 'date', 'ipr_acc', 'ipr_desc', 'go', 'pathways'
        ]
        # Tenta ler. Se tiver mais colunas, o pandas se vira. Se tiver menos, preenche com NaN.
        pfam_df = pd.read_csv(args.pfam, sep='\t', names=col_names, header=None)
        
        # Filtra apenas linhas válidas (remove linhas vazias ou comentários perdidos)
        pfam_df = pfam_df.dropna(subset=['seq_id', 'start_aa', 'end_aa'])
        
    except Exception as e:
        print(f"ERRO ao ler Pfam: {e}")
        return
    print(f"      -> {len(pfam_df)} anotações carregadas.")

    # --- DIAGNÓSTICO DE IDs ---
    gff_ids = set(genes_df['id']).union(set(genes_df['name']))
    pfam_ids = set(pfam_df['seq_id'].astype(str))
    common = gff_ids.intersection(pfam_ids)
    
    if not common:
        print("\n" + "="*60)
        print(" ALERTA: NENHUM MATCH ENTRE GFF E PFAM!")
        print(" O gráfico mostrará apenas os genes (setas), sem domínios.")
        print("-" * 60)
        print(f" IDs no GFF (Exemplo):  {list(gff_ids)[:3]}")
        print(f" IDs no Pfam (Exemplo): {list(pfam_ids)[:3]}")
        print("="*60 + "\n")
        
        # TENTATIVA DE CORREÇÃO AUTOMÁTICA (Se houver apenas 1 gene e 1 proteína)
        if len(gff_ids) == 1 and len(pfam_ids) == 1:
            unica_prot_pfam = list(pfam_ids)[0]
            unico_gene_gff = list(genes_df['id'])[0]
            print(f" -> Detectado caso único. Forçando match: {unica_prot_pfam} assumido como {unico_gene_gff}")
            pfam_df['seq_id'] = unico_gene_gff # Substitui o ID no dataframe do Pfam
    # ---------------------------

    # 3. Mapeamento
    print("[3/4] Construindo gráfico...")
    features = []
    
    # Cores
    color_map = {
        "Pfam": "#ff9999",       # Vermelho claro
        "Phobius": "#99ccff",    # Azul claro
        "Coils": "#cccccc",      # Cinza
        "SignalP": "#ffcc99",    # Laranja
        "TMHMM": "#99ff99"       # Verde
    }
    
    for _, gene in genes_df.iterrows():
        # Adiciona Gene
        features.append(GraphicFeature(
            start=gene['start'], end=gene['end'], strand=gene['strand'],
            color="#e0e0e0", label=gene['name'], thickness=12
        ))
        
        # Busca domínios (Tenta pelo ID e pelo Name)
        doms = pfam_df[pfam_df['seq_id'].astype(str) == str(gene['id'])]
        if doms.empty:
            doms = pfam_df[pfam_df['seq_id'].astype(str) == str(gene['name'])]
            
        for _, d in doms.iterrows():
            # Converte AA -> NT
            # (start_aa - 1) * 3 para 0-based relativo ao inicio do gene
            rel_start = (int(d['start_aa']) - 1) * 3
            rel_end = int(d['end_aa']) * 3
            
            if gene['strand'] == 1:
                g_start = gene['start'] + rel_start
                g_end = gene['start'] + rel_end
            else:
                g_end = gene['end'] - rel_start
                g_start = gene['end'] - rel_end
            
            # Define cor baseada na análise (Pfam, Phobius, etc)
            c = color_map.get(d['analysis'], "#ffd700") # Amarelo default
            
            # Limpa nome do domínio (remove descrição longa se for Phobius/TMHMM)
            label = d['signature_acc'] if pd.notna(d['signature_acc']) else d['analysis']
            if d['analysis'] == 'Pfam': label = d['domain_name']
            
            features.append(GraphicFeature(
                start=g_start, end=g_end, strand=gene['strand'],
                color=c, label=str(label)[:20], # Corta labels muito longos
                thickness=8, fontdict={"fontsize": 6}
            ))

    # 4. Plotar
    record = GraphicRecord(sequence_length=genome_length, features=features)
    ax, _ = record.plot(figure_width=15)
    ax.set_title(args.title, fontsize=16)
    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi)
    print(f"[4/4] Salvo em: {args.output}")

if __name__ == "__main__":
    main()