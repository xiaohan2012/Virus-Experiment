methods=(CEpiMatch MATT Multiprot SPa SPb SPe TMa TMb TMc RMSD)

for method in "${methods[@]}"
do
    m_path="clustering/data/$method.csv"
    
    echo "AUC for $method"
    python yard_gen.py $m_path $method > roc/$method.yard
    yard-auc "roc/$method.yard"
    yard-plot -o roc_plot/$method.png roc/$method.yard
done
