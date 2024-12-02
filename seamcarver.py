#!/usr/bin/env python3

from picture import Picture
import math
import sys
from PIL import Image


class SeamCarver(Picture):
    def energy(self, x: int, y: int) -> float:
        '''
        Return the energy of pixel at column x and row y
        '''
        # Check if the pixel coordinates are valid; raise an error if not
        if not (0 <= x < self.width() and 0 <= y < self.height()):
            raise IndexError("Pixel indices out of range.")

        # Get the pixels directly to the left and right of the current pixel, wrapping around edges.
        left = self[x-1 if x > 0 else self.width() - 1, y]
        right = self[x+1 if x < self.width() - 1 else 0, y]

        # Get the pixels directly above and below the current pixel, wrapping around edges.
        up = self[x, y-1 if y > 0 else self.height() - 1] 
        down = self[x, y+1 if y < self.height() - 1 else 0] 

        # Compute the differences in red, green, and blue components in the x-direction.
        dx_r, dx_g, dx_b = (right[0] - left[0], right[1] - left[1], right[2] - left[2])
        # Compute the differences in red, green, and blue components in the y-direction.
        dy_r, dy_g, dy_b = (down[0] - up[0], down[1] - up[1], down[2] - up[2])

        # This gives us the horizontal and vertical gradient, which measures how much the color changes
        # horizontally and vertically at this pixel. The horizontal and vertical gradients together capture
        # how much the pixel stands out from its surroundings in all directions.

        # Calculate the squared gradients in both x and y directions.
        # These squared gradients represent the magnitude of the color change in each direction, summing
        # across all three color channels (red, green, and blue). The squaring ensures that negative
        # differences (e.g., color becoming darker) don't cancel out positive differences (e.g., color 
        # becoming brighter)
        delta_x2 = dx_r**2 + dx_g**2 + dx_b**2
        delta_y2 = dy_r**2 + dy_g**2 + dy_b**2

        # Compute and return the energy of the pixel as the square root of the sum of squared gradients.
        # The total energy is computed as the square root of the sum of the squared gradients. This step
        # uses the Pythagorean theorem, treating the x and y gradients as perpendicular components of a
        # vector. The result gives the magnitude of the color change at the pixel, regardless of the 
        # direction. This value represents how "important" the pixel is to the structure of the image.
        # High-energy pixels occur at edges or regions with significant color transitions, which should
        # be preserved during seam carving.
        return math.sqrt(delta_x2 + delta_y2)

    def find_vertical_seam(self) -> list[int]:
        '''
        Return a sequence of indices representing the lowest-energy vertical seam
        '''
        width, height = self.width(), self.height()

        # Create a 2D list storing the energy values of all pixels in the image.
        energy_map = [[self.energy(i, j) for i in range(width)] for j in range(height)]

        # Create a 2D list to store cumulative energy values for dynamic programming.
        cumulative_energy = [[0] * width for _ in range(height)]
        # Create a 2D list to store the backtracking path.
        path = [[0] * width for _ in range(height)]

        # Initialize the top row
        # The cumulative energy for the top row of the image is simply the energy values of the pixels
        # in that row because there are no rows above it to contribute to their cumulative energy
        
        # The cumulative energy of each pixel in the rows below will depend on the cumulative energy of
        # the pixels directly above it (or diagonally above it).
        for i in range(width):
            cumulative_energy[0][i] = energy_map[0][i]

        # Populate the energy table
        for j in range(1, height):
            for i in range(width):
                # Get the minimum cumulative energy from the previous row's neighboring pixels.
                min_energy = cumulative_energy[j - 1][i]
                
                # Check the left diagonal pixel if it exists.
                if i > 0:
                    # Compare the cumulative energy of the top-left diagonal pixel (cumulative_energy[j - 1][i - 1])
                    # with the current min_energy and update min_energy to the smaller value
                    min_energy = min(min_energy, cumulative_energy[j - 1][i - 1])
                
                # Check the right diagonal pixel if it exists.
                if i < width - 1:
                    # Repeat of left diagonal except but with the right diagonal
                    min_energy = min(min_energy, cumulative_energy[j - 1][i + 1])

                # Update the cumulative energy for the current pixel.
                cumulative_energy[j][i] = energy_map[j][i] + min_energy

                # Track the path to reconstruct the seam later.
                # Records the column index of the pixel in the previous row (j - 1) that contributed
                # the minimum cumulative energy to the current pixel in row j.
                if min_energy == cumulative_energy[j - 1][i]:
                    # If the min energy comes from the pixel above
                    path[j][i] = i
                elif i > 0 and min_energy == cumulative_energy[j - 1][i - 1]:
                    # If the min energy comes from the pixel to the top-left
                    path[j][i] = i - 1
                else:
                    # If the min energy comes from the pixel to the top-right
                    path[j][i] = i + 1

        # Find the column index of the pixel with the minimum cumulative energy in the last row.
        # cumulative_energy[-1]      -> the last row of the cumulative energy table
        # min(cumulative_energy[-1]) -> finds the smallest cumulative energy value in the last row (optimal vert seam)
        # .index()                   -> retrieves the column index of the pixel with this minimum cumulative energy
        # min_index                  -> start point for backtracking
        min_index = cumulative_energy[-1].index(min(cumulative_energy[-1]))
        # Backtrack to determine the vertical seam.
        seam = [0] * height
        for j in range(height - 1, -1, -1):
            seam[j] = min_index             # Record the column index of the seam pixel.
            min_index = path[j][min_index]  # Move to the previous row's seam pixel.

        return seam

    def remove_vertical_seam(self, seam: list[int]):
        '''
        Remove a vertical seam from the picture
        '''
        # Ensure the seam is valid by checking its length and the image width.
        if len(seam) != self.height():
            raise SeamError("Invalid seam length.")
        if self.width() <= 1: # Ensure the image has at least 2 columns.
            raise SeamError("Cannot remove seam from image with width <= 1.")

        # Remove the pixels in the seam, one row at a time.
        for j in range(self.height()):
            # Ensure the seam is continuous by checking the difference between adjacent entries.
            if j > 0 and abs(seam[j] - seam[j - 1]) > 1:
                raise SeamError("Invalid seam: adjacent entries differ by more than 1.")
            # Shift pixels to the left to fill the gap left by the removed seam pixel.
            for i in range(seam[j], self.width() - 1):
                self[i, j] = self[i + 1, j]
            # Delete the last pixel in the row after shifting.
            del self[self.width() - 1, j]

        self._width -= 1

    def find_horizontal_seam(self) -> list[int]:
        '''
        Return a sequence of indices representing the lowest-energy horizontal seam
        '''
        # Transpose the image to treat horizontal seams as vertical seams.
        transposed_image = self.picture().transpose(Image.Transpose.TRANSPOSE)
        # Create a SeamCarver instance for the transposed image.
        transposed_seam_carver = SeamCarver(transposed_image)
        # Find the vertical seam in the transposed image.
        transposed_seam = transposed_seam_carver.find_vertical_seam()
        # Map the transposed seam back to horizontal seam
        horizontal_seam = transposed_seam
        
        return horizontal_seam

    def remove_horizontal_seam(self, seam: list[int]):
        '''
        Remove a horizontal seam from the picture
        '''
        # Transpose the image to swap width and height
        transposed_image = self.picture().transpose(Image.Transpose.TRANSPOSE)
        # Create a SeamCarver instance for the transposed image
        transposed_seam_carver = SeamCarver(transposed_image)
        # Remove the seam by utilizing the remove_vertical_seam
        transposed_seam_carver.remove_vertical_seam(seam)
        # Transpose the image back to its original dimensions
        updated_image = transposed_seam_carver.picture().transpose(Image.Transpose.TRANSPOSE)

        # Update the current object's pixel data directly
        self._width, self._height = updated_image.size
        self.clear()  # Clear the current pixel data
        for j in range(self._height):
            for i in range(self._width):
                self[i, j] = updated_image.getpixel((i, j))

class SeamError(Exception):
    pass
