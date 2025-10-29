#!/bin/bash

classification_tsv="$1"   # путь к файлу Classification.tsv
pfm_dir="$2"              # путь к директории с *.meme файлами
out_dir="$3"              # директория для вывода результатов

# Заголовок сохраняем отдельно
header=$(head -n 1 "$classification_tsv")

# Создаём выходную директорию, если её нет
mkdir -p "$out_dir"

# Основной выходной файл
output_file="$out_dir/Tomtom_results.tsv"

# очищаем файл и пишем заголовок
echo -e "Query_ID\tTarget_ID\tOptimal_offset\tpvalue\tEvalue\tqvalue\tOverlap\tQuery_consensus\tTarget_consensus\tOrientation" > "$output_file"

# Получаем список уникальных суперклассов (без заголовка)
cut -f3 "$classification_tsv" | tail -n +2 | awk '!seen[$0]++' | while IFS= read -r superclass; do
  # Сохраняем отфильтрованные строки в переменную
  data=$(awk -F'\t' -v s="$superclass" '$3==s' "$classification_tsv")
  # Берём второй столбец (ID) из data в массив
  mapfile -t ids < <(echo "$data" | awk -F'\t' '{print $2}')

  # Двойной перебор всех уникальных пар ID
  for ((i=0; i<${#ids[@]}; i++)); do
      id1=${ids[i]}
      for ((j=i+1; j<${#ids[@]}; j++)); do
          id2=${ids[j]}

          # Пример: проверка наличия файлов
          pfm1="${pfm_dir}/${id1}.meme"
          pfm2="${pfm_dir}/${id2}.meme"
		  # временная папка для пары
		  tmp_dir=$(mktemp -d)
		  tomtom -no-ssc -oc "$tmp_dir" -verbosity 1 -min-overlap 3 -dist pearson -evalue -thresh 20.0 -time 300 "$pfm1" "$pfm2"
		  # объединяем результаты (например, tomtom.txt) в один файл
          if [[ -f "$tmp_dir/tomtom.tsv" ]]; then
            sed -n '2p' "$tmp_dir/tomtom.tsv" >> "$output_file"
          fi
		  
		  # удаляем временную папку
          rm -rf "$tmp_dir"
      done
  done
done