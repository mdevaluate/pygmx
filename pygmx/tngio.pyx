
import numpy as np

cimport numpy as np
#from libc.stdlib cimport free

cdef extern from "tng/tng_io_fwd.h":
    cdef struct tng_trajectory:
        pass

    ctypedef tng_trajectory *tng_trajectory_t


cdef extern from "tng/tng_io.h":
    ctypedef int int64_t
    cdef int64_t TNG_TRAJ_POSITIONS "TNG_TRAJ_POSITIONS"
    cdef int64_t TNG_TRAJ_VELOCITIES "TNG_TRAJ_VELOCITIES"
    cdef int64_t TNG_TRAJ_FORCES "TNG_TRAJ_FORCES"

    ctypedef enum tng_function_status:
         TNG_SUCCESS, TNG_FAILURE, TNG_CRITICAL

    tng_function_status tng_util_trajectory_open(const char *filename, const char mode, tng_trajectory_t *tng_data_p)

    tng_function_status tng_util_trajectory_close(tng_trajectory_t *tng_data_p)

    tng_function_status tng_util_time_of_frame_get(const tng_trajectory_t tng_data, const int64_t frame_nr, double *time)

    tng_function_status tng_util_pos_read_range(const tng_trajectory_t tng_data,
                                                const int64_t first_frame,
                                                const int64_t last_frame,
                                                float **positions,
                                                int64_t *stride_length)

    tng_function_status tng_util_vel_read_range(const tng_trajectory_t tng_data,
                                                const int64_t first_frame,
                                                const int64_t last_frame,
                                                float **velocities,
                                                int64_t *stride_length)

    tng_function_status tng_util_force_read_range(const tng_trajectory_t tng_data,
                                                  const int64_t first_frame,
                                                  const int64_t last_frame,
                                                  float **forces,
                                                  int64_t *stride_length)

    tng_function_status tng_util_box_shape_read_range(const tng_trajectory_t tng_data,
                                                      const int64_t first_frame,
                                                      const int64_t last_frame,
                                                      float **box_shape,
                                                      int64_t *stride_length)

    tng_function_status tng_num_frames_get(const tng_trajectory_t tng_data, int64_t *n)

    tng_function_status tng_num_particles_get(const tng_trajectory_t tng_data, int64_t *n);

    tng_function_status tng_medium_stride_length_get(const tng_trajectory_t tng_data, int64_t *len);

    tng_function_status tng_long_stride_length_get(const tng_trajectory_t tng_data, int64_t *len);

    tng_function_status tng_util_trajectory_next_frame_present_data_blocks_find(
                    const tng_trajectory_t tng_data,
                    int64_t current_frame,
                    const int64_t n_requested_data_block_ids,
                    const int64_t *requested_data_block_ids,
                    int64_t *next_frame,
                    int64_t *n_data_blocks_in_next_frame,
                    int64_t **data_block_ids_in_next_frame)

    tng_function_status tng_util_particle_data_next_frame_read(
                    const tng_trajectory_t tng_data,
                    const int64_t block_id,
                    void **values,
                    char *data_type,
                    int64_t *retrieved_frame_number,
                    double *retrieved_time)


cdef check_status(status, message=''):
    if status != TNG_SUCCESS:
        print(message)

#cdef ptr_to_array(float *ptr, np.npy_intp *shape, int np):
#    arr = np.PyArray_SimpleNewFromData(np, shape, np.NPY_FLOAT, ptr)
#    return arr

cdef class TNGTrajectory:
    cdef:
        tng_trajectory_t trajectory
        int nframes, natoms, stride

    def get_position_interval(self):
        cdef:
            double time
            char data_type
            int frame1, frame2

        cdef np.ndarray[float, ndim=2] positions = np.empty((self.natoms, 3), dtype=np.float32)
        cdef void *data = <void *>positions.data

        tng_util_particle_data_next_frame_read(self.trajectory, TNG_TRAJ_POSITIONS, &data, &data_type, &frame1, &time)
        tng_util_particle_data_next_frame_read(self.trajectory, TNG_TRAJ_POSITIONS, &data, &data_type, &frame2, &time)
        del positions
        return frame2 - frame1

    def __len__(self):
        return self.nframes

    def __getitem__(self, frame):
        cdef np.ndarray[float, ndim=2] positions = np.empty((self.natoms, 3), dtype=np.float32)
        cdef float *pos = <float *>positions.data

        #check_status(
        if (TNG_SUCCESS !=
            tng_util_pos_read_range(self.trajectory, frame, frame, &pos, &self.stride)#,
            #'Reading of frame not successful.'
            ):
            raise IndexError()

        return positions

    def __init__(self, filename):
        assert filename, 'Filename must not be null.'
        if isinstance(filename, str):
            filename = filename.encode()

        if (tng_util_trajectory_open(filename, b'r', &self.trajectory) != TNG_SUCCESS):
            raise OSError('Trajectory could not be opened.')

        tng_num_frames_get(self.trajectory, &self.nframes)
        tng_num_particles_get(self.trajectory, &self.natoms)


from math import ceil

class TNGReader:
    def __init__(self, filename):
        self.trajectory = TNGTrajectory(filename)
        self.nstxout = self.trajectory.get_position_interval()

    def __len__(self):
        return ceil(len(self.trajectory) / self.nstxout)

    def __getitem__(self, item):
        try:
            return self.trajectory[item * self.nstxout]
        except IndexError:
            raise IndexError('Trajectory frame {} out of range.'.format(item))
