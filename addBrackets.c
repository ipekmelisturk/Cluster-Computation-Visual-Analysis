#include <stdio.h>
#include <stdlib.h>

void main()
{
    FILE *file = fopen("output.json", "r+");
    if (file == NULL)
    {
        printf("Failed to open file\n");
    }

    // Move the file position to the end
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);

    // Allocate memory to store the existing content
    char *buffer = (char *)malloc(file_size + 1);
    if (buffer == NULL)
    {
        printf("Failed to allocate memory\n");
        fclose(file);
    }

    // Move the file position to the beginning
    fseek(file, 0, SEEK_SET);

    // Read the existing content into buffer
    fread(buffer, file_size, 1, file);
    buffer[file_size] = '\0';

    // Move the file position to the beginning
    fseek(file, 0, SEEK_SET);

    // Write '[' at the beginning
    fprintf(file, "[");

    // Write the existing content back to the file
    fprintf(file, "%s", buffer);

    // Move the file position to the 3rd last position
    fseek(file, -3, SEEK_END);

    // Write ']' at the 3rd last position
    fprintf(file, "]");

    // Free the memory and close the file
    free(buffer);
    fclose(file);
}