# Example (fill in actual values from literature or assumptions)
p_vision_loss <- 0.01  # Per cycle in Relapse
p_stroke <- 0.005

p_diabetes <- 0.008
p_cvd <- 0.007
p_infection <- 0.01
p_fracture <- 0.006

# Array to track events by state, cycle, and arm
event_array <- function(event_name) {
  array(0, dim = c(n_states, n_cycles, n_arms),
        dimnames = list(state = state_names, cycle = NULL, arms = arms))
}

events <- list(
  vision_loss = event_array("vision_loss"),
  stroke = event_array("stroke"),
  diabetes = event_array("diabetes"),
  cvd = event_array("cvd"),
  infection = event_array("infection"),
  fracture = event_array("fracture")
)

for (a in 1:n_arms) {
  for (c in 2:n_cycles) {
    for (s in 1:n_states) {
      state_name <- state_names[s]
      n_people <- pop[state_name, c - 1, a]
      
      if (state_name != "Dead" && n_people > 0) {
        # GCA complications: only in Relapse
        if (state_name == "Relapse") {
          events$vision_loss[state_name, c, a] <- rbinom(1, n_people, p_vision_loss)
          events$stroke[state_name, c, a]       <- rbinom(1, n_people, p_stroke)
        }
        
        # GC-related AEs: all states except death
        events$diabetes[state_name, c, a]   <- rbinom(1, n_people, p_diabetes)
        events$cvd[state_name, c, a]        <- rbinom(1, n_people, p_cvd)
        events$infection[state_name, c, a]  <- rbinom(1, n_people, p_infection)
        events$fracture[state_name, c, a]   <- rbinom(1, n_people, p_fracture)
      }
    }
  }
}

# Example: cycle-level cost calculation
event_costs <- c(
  vision_loss = 5000,
  stroke = 8000,
  diabetes = 200,
  cvd = 1000,
  infection = 500,
  fracture = 700
)

# Inside your cost loop:
cycle_costs[a, c] <- base_state_cost +
  events$vision_loss[, c, a] %*% rep(event_costs["vision_loss"], n_states) +
  events$stroke[, c, a] %*% rep(event_costs["stroke"], n_states) +
  events$diabetes[, c, a] %*% rep(event_costs["diabetes"], n_states) +
  events$cvd[, c, a] %*% rep(event_costs["cvd"], n_states) +
  events$infection[, c, a] %*% rep(event_costs["infection"], n_states) +
  events$fracture[, c, a] %*% rep(event_costs["fracture"], n_states)
