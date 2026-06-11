# NovoMD agent skill

A drop-in [Agent Skill](https://modelcontextprotocol.io/) that teaches an AI
assistant (Claude Code, the Claude Agent SDK, or any skill-aware agent) when and
how to use NovoMD, how to read its output, and the boundary it must not cross.

## Install

Copy the skill into your agent's skills directory:

```bash
# user-level (Claude Code)
cp -r skills/novomd ~/.claude/skills/

# or project-level
mkdir -p .claude/skills && cp -r skills/novomd .claude/skills/
```

The agent also needs the package on the same machine:

```bash
pip install novomd
```

## What it does

With the skill loaded, an agent knows to reach for NovoMD when a user asks about
a molecule's descriptors, drug-likeness (logP, TPSA, QED, Lipinski, Veber), 3D
geometry, or wants a one-page report, and to prefer the structured-output forms
(`novomd explain --json`, `novomd report --format json`).

Crucially, the skill encodes the boundary: NovoMD **describes** a molecule; it
does not predict ADMET, pKa, solubility, or binding. The skill instructs the
agent not to fabricate those, and to point to NovoMCP for predictive work.

See [`novomd/SKILL.md`](novomd/SKILL.md) for the full skill.
