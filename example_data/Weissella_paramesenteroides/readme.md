This is an example dataset for *Weissella paramesenteroides*.

This directory includes:
- Gene annotation information  
- Gene name modifications  
- Upregulated genes from transcriptomic data  
- The generated metabolic network file: `wp_network.txt`  
- The resulting metabolic gene clusters: `wp_cluster.txt`
- In central.txt, only the central mass values are provided, without retention time.

In the metabolic gene cluster file, an entry like “115.0399-H” represents a negative ion mode, where the actual neutral mass should be calculated as “115.0399 + 1.007825”, with 1.007825 being the exact molecular weight of H⁺ (proton).
