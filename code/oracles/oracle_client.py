#/usr/bin/env python3

import socket
import secrets

from queue import Queue
from threading import Thread
from oracle_helpers import *

class OracleClient():
    def __init__(self, host, port, logging):
        self.port           = port
        self.host           = host
        self.logging        = logging
        self.bufsize        = BUF_SIZE
        self.answer_queue   = Queue()

        self.ans_parsers = {
            OP_ANS_KEYGEN: self._keygen_ansparse,
            OP_ANS_SIGN: self._sign_ansparse,
            OP_ANS_VERIFY: self._verify_ansparse,
            OP_ANS_WRAP: self._wrap_ansparse
        }

    def _receive(self, expected_answers=-1):
        log_tag = f"[OracleClient._receive] "
        self.logging.info(f"{log_tag}Start running client receiver thread")

        data = b""
        while expected_answers != 0:
            try:
                while b"\n" not in data:
                    data += self.socket.recv(self.bufsize)
                
                idx = data.index(b"\n")
                curr_data = data[:idx]
                data = data[idx+1:]

                if not curr_data or len(curr_data) <= 0:
                    break

                parsed = parse_query(curr_data)

                if parsed == None:
                    self.logging.warning(f"Couldn't parse query '{curr_data}', ignoring...")
                    continue
                    
                op, ident, ans_data = parsed

                if ident == OP_CLOSE_CONN:
                    self._close()
                    break

                if op not in self.ans_parsers:
                    self.logging.warning(f"{log_tag}Received unexpected op code {op} in query {curr_data}, when expecting answer, ignoring...")
                    continue

                answer = self.ans_parsers[op](ans_data)
                self.answer_queue.put((ident, answer))
                self.logging.info(f"{log_tag}Received following answer for identifier {ident}: {answer}")
            except Exception as e:
                self.logging.error(f"{log_tag}Caught the following exception, stopping oracle client receiver.\n{e}")
                break

            expected_answers -= 1

        self.logging.info(f"Stopped oracle client receiver")

    def _send_query(self, op, ident, query):
        self.logging.info(f"sending query for operation '{op}', ID '{ident}', and query '{query}'")
        query_enc = build_query(op, ident, query)
        self.socket.sendall(query_enc + b"\n")

    def _connect(self):
        self.logging.info(f"Connect oracle client to {self.host}:{self.port}")
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect((self.host, self.port))

    def _close(self):
        self.socket.close()

    def _send_fin(self):
        self._send_query(OP_CLOSE_CONN, 0, 0)

    def query_all(self, op, queries):
        self.logging.info("Oracle client starting")
        sent_queries = dict()

        try:
            self._connect()

            nqueries = len(queries)
            receiver = Thread(target=self._receive,
                                kwargs={"expected_answers": nqueries})
            receiver.start()

            off = secrets.randbelow(2**32)
            for idx, query in enumerate(queries):
                ident = (off + idx) % (2**32)
                sent_queries[ident] = query
                self._send_query(op, ident, query)

            answer_idents = []
            while len(answer_idents) < nqueries:
                ident, answer = self.answer_queue.get()
                answer_idents.append(ident)
                yield (sent_queries[ident], answer)

            receiver.join()
        except KeyboardInterrupt:
            self.logging.warning("Caught user interrupt, terminating client...")
        finally:
            self._send_fin()
            self._close()

        for ident in sent_queries.keys():
            if ident not in answer_idents:
                self.logging.warning(f"Answer with identifier {ident} missing (for query {sent_queries[ident]}).")
                continue

    def query_single(self, op, query):
        """
        Get the oracle server response for a single query 'query' (must be bytes or int)
        """

        answer = None

        try:
            self._connect()

            ident = secrets.randbelow(2**32)
            self._send_query(op, ident, query)
            self._receive(expected_answers=1)
            answered_ident, answer = self.answer_queue.get()

            if answered_ident != ident:
                self.logging.error(f"Did not receive answer for {ident} instead got answer to: {answered_ident}")
                return None
        except KeyboardInterrupt:
            self.logging.warning("Caught user interrupt, terminating client...")
        finally:
            self._send_fin()
            self._close()

        return answer

    def _keygen_ansparse(self, ans):
        try:
            N = int(ans, 16)
            self.logging.info(f"Received public RSA key modulus N: {N}")
        except Exception as e:
            self.logging.error(f"Failed to parse keygen answer {ans} with error:\n{e}")
            exit(1)

        return N

    def build_keygen_query(self, pub_exp, modulus_bits, label):
        """Query format: pub_exp:modulus_bits:label"""
        return f"{pub_exp}:{modulus_bits}:{label}"

    def keygen(self, query):
        raise NotImplementedError("Must be implemented by child class.")

    def _sign_ansparse(self, ans):
        try:
            sig = int(ans, 16)
            self.logging.info(f"Received RSA signature: {sig}")
        except Exception as e:
            self.logging.error(f"Failed to parse signature {ans} with error:\n{e}")
            exit(1)

        return sig

    def build_sign_query(self, data, label):
        """Query format: hex(data):label"""
        return f"{i2h(data)}:{label}"

    def _verify_ansparse(self, ans):
        try:
            decision = int(ans, 16)
            if decision > 1:
                self.logging.error(f"Received invalid verification decision {decision}")
            decision_bool = (decision == 1)
            self.logging.info(f"Received RSA signature verification result: {decision_bool}")
        except Exception as e:
            self.logging.error(f"Failed to parse signature {ans} with error:\n{e}")
            exit(1)

        return decision_bool

    def build_verify_query(self, data, sig, label):
        """Query format: hex(data):hex(sig):label"""
        return f"{i2h(data)}:{i2h(sig)}:{label}"

    def _wrap_ansparse(self, ans):
        return ans

    def build_wrap_query(self, wrapping_key_label, wrapped_key_label):
        """Query format: wrapping_key_label:wrapped_key_label"""
        return f"{wrapping_key_label}:{wrapped_key_label}"