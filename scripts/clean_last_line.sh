#!/bin/bash

file="$1"

if [ ! -f "$file" ]; then
  echo "❌ 文件不存在：$file"
  exit 1
fi

last_line=$(tail -n 1 "$file")

if [[ -z "$last_line" ]]; then
  echo "🔍 检测到最后一行是空行，正在删除..."
  # 删除最后一行（如果是空的）
  # macOS BSD sed
  if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' '${/^$/d;}' "$file"
  else
    # Linux GNU sed
    sed -i '${/^$/d;}' "$file"
  fi
  echo "✅ 已删除最后的空行"
else
  echo "✅ 最后一行不是空行，无需处理"
fi
