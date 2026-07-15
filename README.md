# reorientation_sequence_impulse_planner

# ChuteTransitionSolver

This program analyzes polyhedron transitions in a two-plane chute setup.
It maps out every pose a rigid part can stably rest in as it moves down
(or is jostled up) a wall-and-floor chute, and lets you:

- map the **optimal transition path** between two specific poses,
- identify every pose **achievable given a transition sequence**
  constraint (e.g. "must alternate Wall → Floor → Wall"), and
- **generate optimal transition sequences** automatically, either for
  one specific pose (**Individual** scope) or across the whole set of
  poses at once (**Global** scope).

## Requirements

- **MATLAB** R2023a or later — the tool uses the `fegeometry` class,
  which was introduced in R2023a.
- **Partial Differential Equation Toolbox** — required for `fegeometry`
  and `generateMesh`, used to import and mesh the part's STL.
- No other toolboxes are required. All graph/geometry logic (convex
  hulls, quaternions, Dijkstra/BFS, sequence search) is implemented from
  scratch on top of base MATLAB.
- *(Optional, non-blocking)* `findjobj` on the path — used only to
  auto-scroll the info panel to the top after an update. It's wrapped in
  `try/catch`, so the tool runs fine without it.
- A part's **`.stl`** file.

## Operation

1. **Select the STL file** for the part you want to analyze.
2. **Enter the chute's angle combination**:
   - **Alpha** = pitch — the chute's tilt along its length (how steep
     the down-chute slide is).
   - **Beta** = roll — the chute's tilt across its width (how much the
     wall leans relative to vertical).
3. The console logs progress as the tool works through its pipeline
   (candidate poses found, stable poses kept, edges built, etc.), and
   the interactive network UI opens automatically once the transition
   graph is complete.

From there you're in the interactive tool — everything below is
available directly in that window.

### The pipeline (what happens before the UI opens)

1. **Mesh & candidate poses** — the STL is meshed, and every face the
   part's convex hull could physically rest on (at any rotation about
   the vertical axis) is enumerated as a candidate pose.
2. **Stability filtering** — a pose survives only if the part's center
   of mass falls inside its support polygon, checked twice: once
   against straight-down gravity, once against the tilted floor/wall
   frame defined by your alpha/beta angles.
3. **Instability scoring** — surviving poses get a moment-arm ratio
   (against both the wall and the floor) that measures how close they
   are to tipping on their own. Poses past a fixed threshold are marked
   as immediate transition sources rather than resting states.
4. **Transition simulation** — from each stable pose, the tool
   physically simulates a tip in each of four directions and follows it
   (chaining through any intermediate unstable poses) until it lands on
   another stable pose:
   - **Down-Wall** / **Down-Floor** — tipping in the direction of chute
     travel, pivoting off the wall or the floor contact edge.
   - **Up-Wall** / **Up-Floor** — tipping against the direction of
     travel (e.g. from vibration, backward handling, or manual reset).
5. **Graph assembly** — every resolved transition becomes a directed,
   costed edge between two stable poses. This graph is what the UI
   queries.

### Network view

The main panel shows every stable pose as a numbered, color-coded
bubble (color reflects its wall/floor stability margin, and whether
it's a dead-end source with no viable further transition). Edges
between bubbles are color-coded by transition type (Down-Wall,
Down-Floor, Up-Wall, Up-Floor), with dashed segments marking "free"
intermediate hops the part passes through on its way to a stable
landing pose.

### Pose inspection

Click any bubble to open a 3D render of that pose — contact points,
centroid, and the part mesh itself — as the **START** pose. Click a
second bubble to render it as the **END** pose. This works alongside
every other mode below; it's how you visually confirm what pose a
number actually represents.

### Path finding (two specific poses)

With a START and END pose selected (and no sequence constraint active),
the tool computes and draws:
- the **fewest-hops** path (breadth-first search), and
- the **lowest-cost** path (Dijkstra, using the cumulative-rotation cost
  model below).

If they're the same path it's drawn once; if they differ, both are
drawn simultaneously so you can compare the trade-off directly.

### Sequence constraints (achievable poses under a rule)

Six dropdown slots let you build a required transition pattern out of
**Wall** and **Floor** categories — e.g. `Wall → Floor → Wall`. Up and
Down variants of the same category are interchangeable within a slot
(an Up-Wall hop can satisfy a "Wall" slot just as well as a Down-Wall
hop). Given that pattern, the tool searches the graph and shows you
every pose that's reachable while honoring it:

- **Individual scope** — search from (or to) one specific anchor pose.
- **Global scope** — run the same search from every pose at once, so
  you see the full reachability picture for that sequence across the
  whole part.
- **Direction** toggle — search `Paths → END` (which poses can reach
  the anchor) or `Paths FROM START` (which poses the anchor can reach).

### Generate Sequences (auto-search)

Rather than hand-building a sequence, this button enumerates every
Wall/Floor category sequence up to length 4, evaluates each one's
coverage (how many poses it connects) and cost/hop efficiency for the
current scope and direction, and loads the best one straight into the
sequence slots — showing its ranked debug output before applying it.

### Filters

- **Allow Types** — restrict which of the four transition types
  (Down-Wall, Down-Floor, Up-Wall, Up-Floor) are eligible for any
  search, path, or sequence query.
- **Max Δ°** — caps the maximum rotation angle allowed for a single
  physical hop. This is applied live at query time, so you can adjust
  it without rebuilding the graph.
- **Metric** — choose whether path/sequence search optimizes for
  **Lowest Cost** or **Fewest Hops**.

### Cost model

- **Down-chute transitions** (Down-Wall, Down-Floor): cost is purely
  the cumulative rotation angle of the physical tip, normalized to
  `[0, 1]`.
- **Up-chute transitions** (Up-Wall, Up-Floor): cost adds a penalty for
  the *source* pose's own instability margin on top of the angle cost —
  tipping against gravity from an already-marginal pose is a real risk,
  not a benefit, so it's weighted accordingly rather than treated the
  same as a down-chute tip.

## Notes

- The transition-type numbering (`1=Down-Wall, 2=Down-Floor, 3=Up-Wall,
  4=Up-Floor`) is a locked convention referenced throughout the code —
  don't renumber it without updating every switch/mod-based lookup.
- `Up-___` and `Down-___` reflects the direction of transition, where
  `Up-___` refers to the part transitioning opposing the sliding direction
  and `Down-___` refers to the part transitioning with the sliding direction
- `__-Wall` and `__-Floor` reflects the type of transition, where
  `__-Wall` refers to a change in wall contact vertices and `__-Floor`
  refers to a change in floor contact vertices
- Large or complex STLs with many resting faces can take a while during
  the candidate-pose sampling pass; console output shows progress at
  each pipeline stage so you can gauge where time is going.
