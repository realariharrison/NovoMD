# Hugging Face Spaces Dockerfile for NovoMD (Gradio + MCP).
# python:3.11-slim + pip (RDKit ships manylinux wheels) — same pattern as the
# repo's root Dockerfile. Replaces the deprecated continuumio/miniconda3 image.
FROM python:3.11-slim

WORKDIR /app

# System libraries RDKit's 2D drawing needs at runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Python deps (gradio[mcp], novomd, and its rdkit/numpy/scipy wheels)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Application
COPY app.py .

# Non-root user for Hugging Face Spaces
RUN useradd -m -u 1000 user
USER user

EXPOSE 7860

ENV GRADIO_SERVER_NAME="0.0.0.0" \
    GRADIO_SERVER_PORT="7860" \
    GRADIO_MCP_SERVER="true"

CMD ["python", "app.py"]
