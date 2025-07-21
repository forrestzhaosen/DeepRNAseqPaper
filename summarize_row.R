# Function to create summary string for a row
summarize_row <- function(row) {
  # Create a string for each column in the format "ColumnName=Value"
  summaries <- sapply(names(row), function(col_name) {
    paste(col_name, "=", row[[col_name]])
  })
  
  # Combine these strings into a single comma-separated string
  paste(summaries, collapse = ", ")
}

# Apply this function to each row of the dataframe
df$Summary = apply(df, 1, summarize_row)

# View the updated dataframe
print(df)