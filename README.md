# StateCharts.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/StateCharts.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/StateCharts.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/StateCharts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/StateCharts.jl)
[![CI/CD](https://github.com/AlgebraicJulia/StateCharts.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/StateCharts.jl/actions/workflows/julia_ci.yml)

A package for working with [state charts](https://en.wikipedia.org/wiki/State_diagram).


## Remaining repo set up tasks
### CodeCov

Set up Codecov credentials for code coverage (If you have trouble, reach out to an AlgebraicJulia organization owner to help with this)

   1. Log into [Codecov](https://codecov.io) with your GitHub account (this requires that you are a member of the AlgebraicJulia organization)
   2. Navigate to the [AlgebraicJulia organization](https://app.codecov.io/gh/AlgebraicJulia)
   3. Select your new repository from the list (e.x. "AlgebraicX")
   4. Note down the `CODECOV_TOKEN` value (It may be in the "Settings" tab if it doesn't show up immediately)
   5. Navigate back to your new GitHub repository and go to the Settings tab
   6. Go to "Security", "Secrets and variables", and "Actions" and click the "New repository secret" button
   7. Give the secret name `CODECOV_TOKEN` and the Secret value is the value you noted from the Codecov settings
   8. Click "Add secret"

### Buildkite

AlgebraicJulia uses [Buildkite](https://buildkite.com/) to submit resource-intensive processes such as building documentation and executing tests to the [HiPerGator](https://www.rc.ufl.edu/about/hipergator/) computing cluster.

While this template comes with a preconfigured `.buildkite/pipeline.yml` file, this repository is not integrated with Buildkite by default. If you would like your repository to use Buildkite to run processes on HiPerGator, tag an issue with @AlgebraicJulia/SysAdmins. 
